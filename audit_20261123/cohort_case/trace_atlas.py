"""Trace numeric expression in a deposited atlas row back to original stable-ID rows."""
from collections import Counter, defaultdict
from itertools import combinations
import csv
import gzip
import json
from pathlib import Path
import tempfile

import numpy as np
from scipy.io import mmread

from cache_entry import cache
from derive_features import feature_records, OUT

TEMP = Path(tempfile.gettempdir())
OBJECTS = ("B", "mast", "ly", "scFAP", "adipo", "endo")
NONZERO = 3  # minimum informative barcodes per source and output row
MIN_CONTROLS = 4


def write_table(name, rows, fields=None):
    with (OUT/name).open("w", encoding="utf8",newline="") as f:
        w = csv.DictWriter(f,fields or list(rows[0]),delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def main():
    with (OUT/"deposit_inventory.tsv").open(encoding="utf8") as f:
        inventory = {r["path"]:r for r in csv.DictReader(f,delimiter="\t")}
    controls = (OUT/"control_symbols.txt").read_text().splitlines()
    panels = feature_records()
    potential = defaultdict(set)
    for rows in dict(panels.values()).values():
        for gene, symbol, _ in rows:
            potential[symbol].add(gene)
    candidates = {s:ids for s,ids in potential.items() if len(ids)>1}
    metadata = {}
    frequency = Counter()
    for sample in panels:
        if not sample.startswith("Massier"):
            continue
        file = cache(inventory[sample+"/barcodes.tsv.gz"])
        with gzip.open(file,"rt") as f:
            cells = [v.strip() for v in f]
        frequency.update(cells)
        metadata[sample]=cells
    objects = {}
    for obj in OBJECTS:
        stem = TEMP/f"hwat_extract_{obj}"
        with Path(f"{stem}_cells.tsv").open(encoding="utf8") as f:
            cells=list(csv.DictReader(f,delimiter="\t"))
        labels = Path(f"{stem}_features.txt").read_text().splitlines()
        values = mmread(f"{stem}_counts.mtx").tocsr()
        objects[obj] = (cells, {s:i for i,s in enumerate(labels)}, values)
    verified, evidence = [], []
    for sample, barcodes in metadata.items():
        print(sample,flush=True)
        index = {b:i for i,b in enumerate(barcodes) if frequency[b]==1}
        matrixfile = cache(inventory[sample+"/matrix.mtx.gz"])
        with gzip.open(matrixfile,"rb") as f:
            counts = mmread(f).tocsr()
        rows = panels[sample][1]
        assert counts.shape == (len(rows),len(barcodes))
        by_id=defaultdict(list)
        by_symbol=defaultdict(list)
        for i,(gene,symbol,_) in enumerate(rows):
            by_id[gene].append(i)
            by_symbol[symbol].append(i)
        if any(len(by_symbol[s])!=1 for s in controls):
            raise ValueError("Non-unique control feature")
        for obj,(cells,indices,values) in objects.items():
            pairs = [(j,index[c["cell"].split("_")[0]]) for j,c in enumerate(cells)
                     if c["cell"].split("_")[0] in index]
            if not pairs:continue
            ci, ri = np.asarray(pairs).T
            used_controls = [s for s in controls if s in indices]
            cv = values[[indices[s] for s in used_controls]][:,ci].toarray()
            rv = counts[[by_symbol[s][0] for s in used_controls]][:,ri].toarray()
            anchor_exact = np.all(np.isclose(rv,cv,atol=1e-7,rtol=1e-7),axis=0)
            # A separate sensitivity trace allows fractional stored counts if taking
            # ceilings exactly recovers all independent integer control counts.
            keep = (rv==np.ceil(np.maximum(cv-1e-9,0))).all(axis=0) & ((rv!=0).sum(axis=0)>=MIN_CONTROLS)
            verified.append(dict(object=obj,source=sample,unique_barcodes=len(ci),anchor_verified=int(keep.sum()),n_controls=len(used_controls),exact_among_verified=int((keep & anchor_exact).sum())))
            ci, ri = ci[keep],ri[keep]
            if len(ci)<NONZERO:continue
            for key,ids in sorted(candidates.items()):
                if key not in indices or key not in by_symbol:continue
                target = values[indices[key],ci].toarray().ravel()
                detected = int(np.count_nonzero(target))
                if detected<NONZERO:continue
                # Compare all raw IDs ever seen with this name, including IDs
                # present under an alternative name in this particular panel.
                hits=[]
                for gene in sorted(ids):
                    for i in by_id[gene]:
                        v=counts[i,ri].toarray().ravel()
                        exact=np.allclose(v,target,rtol=1e-7,atol=1e-7)
                        if exact or np.array_equal(v,np.ceil(np.maximum(target-1e-9,0))):
                            hits.append((gene,rows[i][1],"exact" if exact else "ceil_stored"))
                if len(hits)==1:
                    gene,raw_name,mode=hits[0]
                    evidence.append(dict(object=obj,output_label=key,source=sample,stable_id=gene,
                        match_mode=mode,raw_name=raw_name,raw_name_multiplicity=len(by_symbol[key]),competing_ids="|".join(sorted(ids)),
                        n_cells=len(ci),detected_cells=detected,stored_sum=round(float(target.sum()),6)))
    write_table("barcode_verification.tsv",verified)
    write_table("source_vectors.tsv",evidence)
    by_key=defaultdict(list)
    for e in evidence: by_key[e["object"],e["output_label"]].append(e)
    confirmed=[]
    for (obj,label), entries in sorted(by_key.items()):
        for l,r in combinations(entries,2):
            if l["stable_id"]==r["stable_id"]:continue
            confirmed.append(dict(object=obj,output_label=label,source_a=l["source"],mode_a=l["match_mode"],id_a=l["stable_id"],cells_a=l["n_cells"],detected_a=l["detected_cells"],sum_a=l["stored_sum"],raw_name_a=l["raw_name"],multiplicity_a=l["raw_name_multiplicity"],
                source_b=r["source"],mode_b=r["match_mode"],id_b=r["stable_id"],cells_b=r["n_cells"],detected_b=r["detected_cells"],sum_b=r["stored_sum"],raw_name_b=r["raw_name"],multiplicity_b=r["raw_name_multiplicity"]))
    fields="object output_label source_a mode_a id_a cells_a detected_a sum_a raw_name_a multiplicity_a source_b mode_b id_b cells_b detected_b sum_b raw_name_b multiplicity_b".split()
    write_table("confirmed_object_rows.tsv",confirmed,fields)
    print(len(confirmed),"source pairs",len({r["output_label"] for r in confirmed}),"distinct labels",flush=True)

if __name__=="__main__":
    main()

