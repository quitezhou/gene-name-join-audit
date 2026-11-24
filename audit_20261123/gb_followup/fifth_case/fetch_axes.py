"""Pin deposited inputs and extract study feature axes, retaining per-gene spatial UMI."""
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
import csv,gzip,hashlib,io,json,sys,time,threading
from pathlib import Path
import h5py
import numpy as np
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

HERE = Path(__file__).resolve().parent
SRC = HERE/"source"
AXIS_LOCK = threading.Lock()
URL = "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/{n}/E-MTAB-13{n}/Files/"
def endpoint(n, name):
    return URL.format(n=n) + name
def sess():
    s=requests.Session()
    s.headers.update({"Accept-Encoding":"identity"})
    s.mount("https://",HTTPAdapter(max_retries=Retry(total=3,backoff_factor=.5,status_forcelist=[429,500,502,503,504])))
    return s
def digest(b):
    return hashlib.sha256(b).hexdigest()
def items(obj):
    if isinstance(obj,dict):
        if "path" in obj: yield (obj["path"],obj["size"])
        for child in obj.values():yield from items(child)
    elif isinstance(obj,list):
        for child in obj:yield from items(child)
def inventory(n):
    return dict(items(json.loads((SRC/f"E-MTAB-13{n}.json").read_text(encoding="utf-8"))))
def dump_table(file,cols,rows):
    with file.open("w",encoding="utf-8",newline="") as f:
        w=csv.writer(f,delimiter="\t",lineterminator="\n");w.writerow(cols);w.writerows(rows)
def get_all(url,size):
    for trial in range(6):
        try:
            with sess() as s:
                r=s.get(url,headers={"Range":f"bytes=0-{size-1}"},timeout=(15,75))
                r.raise_for_status()
                if len(r.content)==size:return r.content
        except requests.RequestException as exc:
            print("retrieval retry",url.rsplit("/",1)[1],trial+1,type(exc).__name__,flush=True)
        time.sleep(min(15,1+trial*2))
    raise OSError(f"File size differed from manifest for {url}")

def regular(spec):
    n,name,size=spec
    memo=SRC/(name+".manifest.json")
    if memo.exists():return tuple(json.loads(memo.read_text("utf8")))
    b=get_all(endpoint(n,name),size)
    if n=="526":
        axis=gzip.decompress(b)
        if not axis.endswith(b"\n"):axis+=b"\n"
        count=axis.count(b"\n")
        spots=umis=0
    else:
        with h5py.File(io.BytesIO(b),"r") as h:
            g=h["matrix"];ft=g["features"]
            ids=ft["id"].asstr()[:];names=ft["name"].asstr()[:];types=ft["feature_type"].asstr()[:]
            shape=g["shape"][:]
            ix=g["indices"][:];v=g["data"][:];ptr=g["indptr"][:]
            assert np.all(types=="Gene Expression") and len(ids)==shape[0] and len(set(ids))==len(ids)
            assert np.all(v>0)
            totals=np.bincount(ix,weights=v,minlength=len(ids)).astype(np.int64)
            detected=np.bincount(ix,minlength=len(ids))
            assert np.sum(v)==sum(totals) and len(ptr)==shape[1]+1
            axis="".join(f"{i}\t{s}\t{t}\n" for i,s,t in zip(ids,names,types)).encode()
            count,spots=map(int,shape);umis=int(sum(totals))
            np.savez_compressed(SRC/(name+".npz"),total=totals,detected=detected)
    key=digest(axis)
    with AXIS_LOCK:
        axis_file=SRC/(key+".tsv.gz")
        if not axis_file.exists():
            axis_file.write_bytes(gzip.compress(axis,mtime=0))
    record=(n,name,size,digest(b),key,count,spots,umis)
    memo.write_text(json.dumps(record),encoding="utf8")
    print("done",name,flush=True)
    return record
class Remote(io.RawIOBase):
    def __init__(self,url,size,chunk=256*1024):
        self.s=sess();self.url=url;self.size=size;self.chunk=chunk;self.pos=0;self.cache={}
        h=self.s.head(self.url,timeout=25); h.raise_for_status()
        self.etag=h.headers["ETag"]
        assert int(h.headers["Content-Length"])==size
    def seek(self,offset,whence=0):
        self.pos={0:0,1:self.pos,2:self.size}[whence]+offset
        return self.pos
    def tell(self):
        return self.pos
    def readinto(self,b):
        end=min(self.pos+len(b),self.size);n=0
        while self.pos<end:
            k,off=divmod(self.pos,self.chunk)
            if k not in self.cache:
                lo=k*self.chunk;hi=min(lo+self.chunk,self.size)-1
                for t in range(6):
                    try:
                        r=self.s.get(self.url,headers={"Range":f"bytes={lo}-{hi}","If-Match":self.etag},timeout=(15,30))
                        r.raise_for_status()
                    except requests.RequestException:
                        if t==5:raise
                        continue
                    if r.status_code==206 and r.headers.get("Content-Range")==f"bytes {lo}-{hi}/{self.size}" and r.headers.get("ETag")==self.etag and len(r.content)==hi-lo+1:break
                    print("retry range",k,r.status_code,r.headers.get("Content-Range"),len(r.content),flush=True)
                else:raise OSError(f"Range integrity failed {self.url}:{lo}-{hi}")
                self.cache[k]=r.content
            chunk=self.cache[k][off:off+end-self.pos]
            b[n:n+len(chunk)]=chunk;self.pos+=len(chunk);n+=len(chunk)
        return n
    def close(self):
        self.s.close();super().close()
def vector(group,key):
    d=group[key]
    if "__categories" in group and key in group["__categories"]:
        cats=group["__categories"][key].asstr()[:].tolist()
        return [cats[i] if i>=0 else "" for i in d[:]]
    return d.asstr()[:].tolist()
def process_ref(spec):
    name,size=spec
    with Remote(endpoint("526",name),size) as rr,h5py.File(rr,"r") as h:
        meta={"filename":name,"url":rr.url,"etag":rr.etag,"bytes":size,"root":list(h)}
        for seg in ["var","raw/var"]:
            if seg not in h:continue
            g=h[seg];keys=list(g); ix=g.attrs["_index"]
            print(name,seg,keys,flush=True)
            columns=[ix]+[k for k in keys if k.lower() in ["ensembl","gene_ids"]]
            rows=list(zip(*(vector(g,k) for k in columns)))
            out=SRC/(name+"."+seg.replace("/","-")+".tsv")
            dump_table(out,columns,rows);meta[seg]={"rows":len(rows),"columns":keys,"table_sha256":digest(out.read_bytes())}
        obs=h["obs"];meta["cells"]=len(obs[obs.attrs["_index"]])
        if "environment" in obs:meta["environment_counts"]=dict(Counter(vector(obs,"environment")))
        meta["downloaded_blocks"]=len(rr.cache)
    (SRC/(name+".json")).write_text(json.dumps(meta,indent=2),encoding="utf-8")
    print(name,meta,flush=True)
def main():
    if sys.argv[1:] == ["references"]:
        with ThreadPoolExecutor(2) as pool:
            list(pool.map(process_ref,[(k,n) for k,n in inventory("526").items() if k.endswith(".h5ad")]))
    else:
        tasks=[(n,k,v) for n,suffix in [("530",".h5"),("526","features.tsv.gz")] for k,v in inventory(n).items() if k.endswith(suffix)]
        with ThreadPoolExecutor(8) as pool:
            records=list(pool.map(regular,tasks))
        dump_table(HERE/"input_axes.tsv",["accession_tail","file","bytes","sha256","axis_id","features","spots","total_UMI"], sorted(records))

if __name__=="__main__":main()
