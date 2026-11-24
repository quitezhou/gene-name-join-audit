"""Recover stable IDs for stored S42 expression using all spot counts."""
from collections import defaultdict
import csv
from pathlib import Path
from tempfile import gettempdir

import h5py
import numpy as np
from scipy.io import mmread
from scipy.sparse import csc_matrix
from cache_entry import cache
from derive_features import OUT, feature_records
from trace_atlas import write_table

def main():
    with (OUT/"deposit_inventory.tsv").open() as f:
        rec = next(r for r in csv.DictReader(f,delimiter="\t") if r["path"]=="S42_trimmed/filtered_feature_bc_matrix.h5")
    with h5py.File(cache(rec)) as f:
        g=f["matrix"]
        raw=csc_matrix((g["data"][:],g["indices"][:],g["indptr"][:]),shape=g["shape"][:]).tocsr()
        barcodes={b.decode():i for i,b in enumerate(g["barcodes"][:])}
        ids = list(g["features/id"].asstr()[:])
        symbols = list(g["features/name"].asstr()[:])
    assert len(ids)==len(set(ids))
    d={gene:i for i,gene in enumerate(ids)}
    by_label=defaultdict(list)
    for i,s in enumerate(symbols):by_label[s].append(i)
    all_ids=defaultdict(set)
    for rows in dict(feature_records().values()).values():
        for gene,s,_ in rows:all_ids[s].add(gene)
    prefix=Path(gettempdir())/"hwat_extract_spatial_S42"
    names=Path(f"{prefix}_features.txt").read_text().splitlines()
    cells=Path(f"{prefix}_cells.txt").read_text().splitlines()
    idx=[barcodes[c.split("_")[0]] for c in cells]
    assert len(set(idx))==len(idx)
    stored=mmread(f"{prefix}_counts.mtx").tocsr()
    controls=set((OUT/"control_symbols.txt").read_text().splitlines()) & set(names)
    for s in controls:
        target=stored[names.index(s)].toarray().ravel()
        assert np.array_equal(target,raw[by_label[s][0],idx].toarray().ravel()),s
    evidence=[]
    for i,s in enumerate(names):
        if len(all_ids[s])<2:continue
        target=stored[i].toarray().ravel()
        n=int(np.count_nonzero(target))
        if n<3:continue
        hits=[g for g in sorted(all_ids[s]) if g in d and np.array_equal(target,raw[d[g],idx].toarray().ravel())]
        if len(hits)==1:
            evidence.append(dict(object="spatial_baseline_S42",output_label=s,stable_id=hits[0],raw_name=symbols[d[hits[0]]],raw_name_multiplicity=len(by_label[s]),n_spots=len(cells),detected_spots=n,stored_sum=int(target.sum()),controls_verified=len(controls)))
    write_table("spatial_S42_vectors.tsv",evidence)
    print(len(evidence),len(cells),len(controls))
if __name__=="__main__":
    main()
