"""Retrieve only feature axes and group metadata from authors' portal H5ADs using HTTP Range."""
import csv, gzip, io, json, time, hashlib
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import h5py
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

OUT = Path(__file__).resolve().parent / "source"
BASE = "https://cellgeni.cog.sanger.ac.uk/spinal-development/"
NAMES = ["fetal_spine_cellxgeneNov2023", "cervical_and_lumbar_week_7", "cervical_and_thoracic_week_7","cervical_week_9","lumbar_and_thoracic_week_7","lumbar_week_8","lumbar_week_9","thoracic_week_9"]

class RangeReader(io.RawIOBase):
    size_chunk = 1024**2
    def __init__(self, url):
        self.url = url
        self.session = requests.Session()
        self.session.mount("https://", HTTPAdapter(max_retries=Retry(total=4, backoff_factor=0.5, status_forcelist=(429,500,502,503,504))))
        h = self.session.head(url, timeout=30)
        h.raise_for_status()
        self.size = int(h.headers["Content-Length"])
        self.etag = h.headers["ETag"]
        self.pos = 0
        self.cache = {}
    def seek(self, offset, whence=0):
        self.pos = offset + {0:0,1:self.pos,2:self.size}[whence]
        return self.pos
    def tell(self):
        return self.pos
    def readinto(self,b):
        n = min(len(b),self.size-self.pos)
        count = 0
        while count<n:
            idx,offset = divmod(self.pos,self.size_chunk)
            if idx not in self.cache:
                low = idx*self.size_chunk
                high = min(low+self.size_chunk,self.size)-1
                r=self.session.get(self.url,headers={"Range":f"bytes={low}-{high}", "If-Match":self.etag},timeout=65)
                r.raise_for_status()
                assert r.status_code==206 and r.headers["Content-Range"]==f"bytes {low}-{high}/{self.size}" and r.headers["ETag"]==self.etag
                assert len(r.content)==high-low+1
                self.cache[idx]=r.content
            piece = self.cache[idx][offset:offset+n-count]
            b[count:count+len(piece)] = piece
            count += len(piece)
            self.pos += len(piece)
        return n
    def close(self):
        self.session.close()
        super().close()

def fetch(name):
    url = BASE+name+".h5ad"
    for trial in range(3):
        try:
            with RangeReader(url) as src, h5py.File(src, "r") as f:
                v=f["var"]
                names=v[v.attrs["_index"]].asstr()[:].tolist()
                ids=v["gene_ids"].asstr()[:].tolist() if "gene_ids" in v else [""]*len(names)
                obs=f["obs"]
                meta={"url":url,"etag":src.etag,"h5ad_bytes":src.size,"n_obs":len(obs[obs.attrs["_index"]]),"n_vars":len(names),
                      "groups":list(f),"obs_fields":list(obs),"var_fields":list(v),"ordered_groups":{k:list(f[k]) for k in ("uns","obsm","varm","layers") if k in f},
                      "id_symbol_text_sha256":hashlib.sha256("".join(i+"\t"+s+"\n" for i,s in zip(ids,names)).encode()).hexdigest()}
            with gzip.open(OUT/(name+".tsv.gz"),"wt",encoding="utf-8",newline="") as fo:
                w=csv.writer(fo,delimiter="\t",lineterminator="\n");w.writerow(["stable_id","index"]);w.writerows(zip(ids,names))
            (OUT/(name+".json")).write_text(json.dumps(meta,indent=2)+"\n",encoding="utf-8")
            print(name,meta["n_vars"],meta["obs_fields"][:10],meta["id_symbol_text_sha256"],flush=True)
            return
        except Exception as e:
            if trial==2:raise
            print(name,type(e).__name__,str(e)[:130],flush=True)
            time.sleep(1+trial)

if __name__=="__main__":
    OUT.mkdir(exist_ok=True)
    with ThreadPoolExecutor(4) as pool:
        list(pool.map(fetch,NAMES))
