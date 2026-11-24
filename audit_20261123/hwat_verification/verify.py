"""Independently test the hWAT candidate traces against every gene row.

Reads the published deposits and candidate/count exports produced by the sibling
cohort_case scripts. Barcode uniqueness is assessed among the nine Massier input
matrices; the stored orig.ident is retained as an independent descriptive field.
"""
import csv
import gzip
import hashlib
import json
import re
import sys
import tempfile
from collections import Counter
from pathlib import Path

import h5py
import numpy as np
from scipy.io import mmread
from scipy.sparse import csc_matrix

HERE=Path(__file__).resolve().parent
CASE=HERE.parent/"cohort_case"
sys.path.insert(0,str(CASE))
from cache_entry import cache
TMP=Path(tempfile.gettempdir())

def tab(file):
    with file.open(encoding="utf-8",newline="") as f:
        return list(csv.DictReader(f,delimiter="\t"))

def write(name,rows):
    with (HERE/name).open("w",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t",lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

def fingerprint(matrix, target, names, expected, rng):
    # CSR has sorted/unique sparse column indices. Restrict first by equal
    # support size, then verify each index, then test numeric values.
    on=np.flatnonzero(target>0)
    ind=np.where(np.diff(matrix.indptr)==len(on))[0]
    support=[]
    exact=[]
    rounded=[]
    ceil=np.ceil(np.maximum(target-1e-9,0))
    for i in ind:
        b,e=matrix.indptr[i:i+2]
        if not np.array_equal(matrix.indices[b:e],on):continue
        support.append(i)
        v=matrix.data[b:e]
        if np.array_equal(v,target[on]):exact.append(i)
        if np.array_equal(v,ceil[on]):rounded.append(i)
    # Jaccard tests all raw rows even when pre-processing removes detections.
    nnz=np.diff(matrix.indptr)
    overlap=matrix[:,on].getnnz(axis=1)
    j=overlap/(nnz+len(on)-overlap)
    rank=int(np.count_nonzero(j>j[expected])+1)
    tied=int(np.count_nonzero(j==j[expected]))
    other=float(max(np.delete(j,expected)))
    shuffle=rng.permutation(target.size)[:len(on)]
    permoverlap=matrix[:,shuffle].getnnz(axis=1)
    return {"raw_features":matrix.shape[0],"matched_barcodes":target.size,"target_n_detected":len(on),
            "target_exact_informative_entries":int(np.sum((target[on]==matrix[expected,on].toarray().ravel()))),
            "n_exact_rows":len(exact),"exact_id":names[exact[0]] if len(exact)==1 else "",
            "n_ceil_rows":len(rounded),"ceil_id":names[rounded[0]] if len(rounded)==1 else "",
            "n_support_rows":len(support),"support_id":names[support[0]] if len(support)==1 else "",
            "expected_support_jaccard":round(float(j[expected]),6),"best_other_support_jaccard":round(other,6),
            "expected_support_rank":rank,"expected_rank_ties":tied,
            "n_shuffled_perfect_support_rows":int(np.count_nonzero((permoverlap==len(on))&(nnz==len(on))))}

def check():
    inv={r["path"]:r for r in tab(CASE/"deposit_inventory.tsv")}
    expected_src=tab(CASE/"source_vectors.tsv")
    spatial={r['output_label']:r['stable_id'] for r in tab(CASE/"spatial_S42_vectors.tsv")}
    switches={r['output_label'] for r in tab(CASE/"confirmed_object_rows.tsv")}
    cross={r['output_label'] for r in expected_src if r['output_label'] in spatial and spatial[r['output_label']]!=r['stable_id']}
    labels=switches|cross
    target_rows={(r['source'],r['object'],r['output_label']):r['stable_id'] for r in expected_src if r['output_label'] in labels}
    controls=(CASE/"control_symbols.txt").read_text().splitlines()
    rng=np.random.default_rng(20261123)
    bars={}
    for path,entry in inv.items():
        if path.startswith('Massier') and path.endswith('/barcodes.tsv.gz'):
            with gzip.open(cache(entry),'rt') as f:bars[path.rsplit('/',1)[0]]=[line.strip() for line in f]
    frequency=Counter(b for bs in bars.values() for b in bs)
    objects={}
    for obj in sorted({o for s,o,k in target_rows}):
        prefix=TMP/f"hwat_extract_{obj}"
        keys=Path(f"{prefix}_features.txt").read_text().splitlines()
        assert len(set(keys))==len(keys)
        objects[obj]=(tab(Path(f"{prefix}_cells.tsv")),keys,mmread(f"{prefix}_counts.mtx").tocsr())
    traces=[]
    for sample in sorted({s for s,o,l in target_rows}):
        bs=bars[sample]
        by_barcode={b:i for i,b in enumerate(bs) if frequency[b]==1}
        feature_entry=inv[sample+'/features.tsv.gz']
        raw_bytes=(CASE/'source'/(feature_entry['file_id']+'.tsv.gz')).read_bytes()
        assert hashlib.sha256(raw_bytes).hexdigest()==feature_entry['sha256']
        feature=[l.split('\t') for l in gzip.decompress(raw_bytes).decode().splitlines()]
        ids=[g[0] for g in feature]
        names=[g[1] for g in feature]
        assert len(ids)==len(set(ids))
        with gzip.open(cache(inv[sample+'/matrix.mtx.gz']),'rb') as f:src=mmread(f).tocsc()
        assert src.shape==(len(ids),len(bs))
        for obj,(cells,keys,stored) in objects.items():
            pairs=[(j,by_barcode[c['cell'].split('_')[0]]) for j,c in enumerate(cells) if c['cell'].split('_')[0] in by_barcode]
            if not pairs:continue
            oi,si=np.asarray(pairs).T
            # The candidate rows never contribute to the anchor.
            ck=[s for s in controls if s in keys]
            assert set(ck).isdisjoint(labels) and all(names.count(s)==1 for s in ck)
            sv=stored[[keys.index(s) for s in ck]][:,oi].toarray()
            rv=src[[names.index(s) for s in ck]][:,si].toarray()
            number = re.search(r'#(\d+)',sample).group(1)
            label_consistent=np.array([bool(re.search(r'#'+number+r'\b',cells[j]['cohort'])) for j in oi])
            anchors={'barcode_and_metadata':label_consistent, 'barcode_only':np.ones(len(oi),dtype=bool),
                     'control_detection':((rv>0)==(sv>0)).all(0)&((rv>0).sum(0)>=4),
                     'control_ceiling':(rv==np.ceil(np.maximum(sv-1e-9,0))).all(0)&((rv>0).sum(0)>=4)}
            for arm,keep in anchors.items():
                if keep.sum()<3:continue
                m=src[:,si[keep]].tocsr()
                m.sum_duplicates();m.eliminate_zeros();m.sort_indices()
                if np.any(m.data<0):raise ValueError("negative raw measurement")
                for lab in sorted(labels):
                    aid=target_rows.get((sample,obj,lab))
                    if not aid:continue
                    target=stored[keys.index(lab),oi[keep]].toarray().ravel()
                    if (target>0).sum()<3:continue
                    same_symbol=[i for i,n in enumerate(names) if n==lab]
                    others=[i for i in same_symbol if ids[i]!=aid]
                    summed=np.asarray(m[same_symbol].sum(axis=0)).ravel()
                    alternate_detection_cells=int(np.count_nonzero(np.asarray(m[others].sum(axis=0)).ravel())) if others else 0
                    traces.append(dict(raw_same_symbol_count=len(same_symbol), other_same_symbol_detected_cells=alternate_detection_cells,
                                       pooled_raw_support_matches=int(np.array_equal(summed>0,target>0)),
                                       pooled_raw_ceil_matches=int(np.array_equal(summed,np.ceil(np.maximum(target-1e-9,0)))),source=sample,object=obj,label=lab,expected_id=aid,anchor=arm,n_controls=len(ck),
                                       metadata_cohorts=json.dumps(dict(Counter(cells[j]['cohort'] for j in oi[keep])),sort_keys=True),
                                       **fingerprint(m,target,ids,ids.index(aid),rng)))
        print(sample,flush=True)
    write("atlas_all_rows.tsv",traces)

    with h5py.File(cache(inv["S42_trimmed/filtered_feature_bc_matrix.h5"])) as f:
        g=f['matrix']
        m=csc_matrix((g['data'][:],g['indices'][:],g['indptr'][:]),shape=g['shape'][:])
        ids=g["features/id"].asstr()[:].tolist()
        bs={b:i for i,b in enumerate(g['barcodes'].asstr()[:])}
    fixed=TMP/"hwat_extract_spatial_S42"
    keys=Path(f"{fixed}_features.txt").read_text().splitlines()
    cells=Path(f"{fixed}_cells.txt").read_text().splitlines()
    si=[bs[c.split('_')[0]] for c in cells]
    assert len(set(si))==len(si) and len(ids)==len(set(ids))
    m=m[:,si].tocsr()
    m.sum_duplicates();m.eliminate_zeros();m.sort_indices()
    stored=mmread(f"{fixed}_counts.mtx").tocsr()
    sp_rows=[]
    for lab in sorted(cross):
        target=stored[keys.index(lab)].toarray().ravel()
        aid=spatial[lab]
        sp_rows.append(dict(source="S42",label=lab,expected_id=aid,**fingerprint(m,target,ids,ids.index(aid),rng)))
    write("spatial_all_rows.tsv",sp_rows)

if __name__=='__main__':
    check()


