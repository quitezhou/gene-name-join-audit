"""Trace base and dot-suffixed rows directly from the deposited adipocyte RDS."""
from collections import Counter
import gzip
import numpy as np
from scipy.io import mmread
from verify import CASE,HERE,tab,write,cache

def main():
    inv={r['path']:r for r in tab(CASE/'deposit_inventory.tsv')}
    bars={}
    for path,entry in inv.items():
        if path.startswith('Massier') and path.endswith('/barcodes.tsv.gz'):
            with gzip.open(cache(entry),'rt') as f:bars[path.rsplit('/',1)[0]]=[line.strip() for line in f]
    frequency=Counter(b for bs in bars.values() for b in bs)
    source=HERE/'source'
    keys=(source/'adipo_names.txt').read_text().splitlines()
    cells=tab(source/'adipo_cells.tsv')
    stored=mmread(source/'adipo_selected.mtx').tocsr()
    out=[]
    for sample,num in [('Massier et al. (#1)','1'),('Massier et al. (#4)/21L005051','4')]:
        index={b:i for i,b in enumerate(bars[sample]) if frequency[b]==1}
        pairs=[(j,index[c['cell'].split('_')[0]]) for j,c in enumerate(cells) if c['cohort']=='Massier et al. #'+num and c['cell'].split('_')[0] in index]
        orig,raw=np.asarray(pairs).T
        rowfile=inv[sample+'/features.tsv.gz']
        with gzip.open(CASE/'source'/(rowfile['file_id']+'.tsv.gz'),'rt',encoding='utf8') as f:
            rows=[l.rstrip('\r\n').split('\t') for l in f]
        with gzip.open(cache(inv[sample+'/matrix.mtx.gz']),'rb') as f:
            m=mmread(f).tocsc()[:,raw].tocsr()
        m.sum_duplicates();m.eliminate_zeros();m.sort_indices()
        allnnz=m.getnnz(axis=1)
        for j,k in enumerate(keys):
            target=stored[j,orig].toarray().ravel()
            on=np.flatnonzero(target>0)
            hits=[]
            if len(on):
                for i in np.flatnonzero(allnnz==len(on)):
                    start,end=m.indptr[i:i+2]
                    if np.array_equal(m.indices[start:end],on):hits.append(i)
            base=k[:-2] if k.endswith('.1') else k
            ids_with_symbol=[r[0] for r in rows if r[1]==base]
            out.append(dict(sample=sample,stored_k=k,n_cells=len(orig),n_detected=len(on),raw_feature_candidates=len(rows),
                            support_matches=len(hits) if len(on) else 'not_tested_empty',
                            unique_id=rows[hits[0]][0] if len(hits)==1 else '',
                            unique_raw_name=rows[hits[0]][1] if len(hits)==1 else '',source_raw_ids_in_input_order='|'.join(ids_with_symbol)))
    write("base_and_suffix_trace.tsv",out)

if __name__=='__main__':main()
