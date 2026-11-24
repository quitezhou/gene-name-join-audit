"""Independently identify reference rows from per-cell raw counts, across all raw genes."""
import csv,gzip,hashlib,json
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from fetch_axes import HERE,SRC
from fetch_validation_inputs import local

def verify(env):
    w=np.load(SRC/(env+'_ref_window.npz'));sample=str(w["sample"])
    barcodes=[c.split('-')[0]+'-1' for c in w["cell"]]
    assert len(set(barcodes))==len(barcodes)
    wanted={b:i for i,b in enumerate(barcodes)}
    pathb=local(sample+"-barcodes.tsv.gz"); pathm=local(sample+"-matrix.mtx.gz")
    cols={}
    with gzip.open(pathb,"rt") as f:
        for i,line in enumerate(f,1):
            if (b:=line.strip()) in wanted:cols[i]=wanted[b]
    assert set(cols.values())==set(range(len(barcodes))) and len(cols)==len(barcodes)
    with gzip.open(pathm,"rt") as f:
        while (l:=next(f)).startswith('%'):pass
        n_genes, n_barcodes, nnz=map(int,l.split())
        raws=[];samples=[];data=[]
        for block in pd.read_csv(f,sep=r'\s+',names=['g','b','n'],dtype='int32',chunksize=1000000):
            sel=block[block.b.isin(cols)]
            if len(sel):
                raws.append(sel.g.to_numpy()-1); samples.append(sel.b.map(cols).to_numpy());data.append(sel.n.to_numpy())
    m=csr_matrix((np.concatenate(data),(np.concatenate(raws),np.concatenate(samples))),shape=(n_genes,len(barcodes)))
    m.sum_duplicates();m.sort_indices();m.eliminate_zeros()
    assert (m.data>0).all()
    manifest=pd.read_csv(HERE/"input_axes.tsv",sep="\t").set_index('file');key=manifest.loc[sample+"-features.tsv.gz","axis_id"]
    axis=pd.read_csv(SRC/(key+".tsv.gz"),sep="\t",names=["id","name","type"],keep_default_na=False)
    indexed_names=pd.Index(axis.name)
    from anndata.utils import make_index_unique
    lookup={k:i for i,k in enumerate(make_index_unique(indexed_names))}
    out=[]
    for j,s in enumerate(w["labels"]):
        t=w["counts"][:,j];on=np.flatnonzero(t)
        e=lookup[s];match=[]
        if len(on):
            for i in np.flatnonzero(np.diff(m.indptr)==len(on)):
                a,b=m.indptr[i:i+2]
                if np.array_equal(m.indices[a:b],on) and np.array_equal(m.data[a:b],t[on]):
                    match.append(axis.id[i])
        out.append(dict(environment=env,sample=sample,input_label=s,source_axis_id=axis.id[e],cells=len(t),positive=len(on),
                        total_UMI=int(t.sum()),source_numeric_equal=np.array_equal(m[e,:].toarray().ravel(),t),
                        all_row_unique=len(match)==1,exact_candidates=len(match) if len(on) else None,
                        matched_id=match[0] if len(match)==1 else "",raw_genes=axis.shape[0]))
    meta=dict(source_files=[{"file":p.name,"sha256":hashlib.sha256(p.read_bytes()).hexdigest()} for p in [pathb,pathm]],
              window_start=int(w['window_start']),window_size=int(w['window_size']),barcodes=len(barcodes),sample=sample)
    (SRC/(env+'_numeric_check.json')).write_text(json.dumps(meta,indent=2),encoding='utf8')
    return out
if __name__=="__main__":
    rows=[]
    for env in ["Tumour","Background"]:
        rows.extend(verify(env))
    df=pd.DataFrame(rows);df.to_csv(HERE/"reference_count_fingerprints.tsv",sep="\t",index=False)
    print(df.to_string(index=False))
