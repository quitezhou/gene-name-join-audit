"""Audit deposited targeted and parent axes, with an optional nonzero-row sensitivity."""
import argparse
import csv
import gzip
import hashlib
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "systematic_mechanism"))
from audit_join import Axis, analyze, EVENT_FIELDS


def records(path):
    with path.open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def save(path, rows, fields):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def extract(cache):
    import h5py
    import numpy as np
    from scipy.sparse import csc_matrix

    barcodes = {}
    metadata = []
    with gzip.open(HERE / "source_profiles.tsv.gz", "wt", encoding="utf-8", newline="") as f:
        out = csv.writer(f, delimiter="\t", lineterminator="\n")
        out.writerow(["sample","stable_id","symbol","umi","detected_spots"])
        for rec in records(HERE / "inputs.tsv"):
            name = rec["sample"]
            path = cache / (name + ".h5")
            assert hashlib.sha256(path.read_bytes()).hexdigest() == rec["sha256"]
            with h5py.File(path) as h5:
                g = h5["matrix"]
                ids, symbols = [g[f"features/{k}"].asstr()[:] for k in ("id","name")]
                barcodes[name] = g["barcodes"].asstr()[:].tolist()
                assert len(set(barcodes[name])) == len(barcodes[name])
                assert set(g["features/feature_type"].asstr()[:]) == {"Gene Expression"}
                features = g["features"]
                panels = list(features["target_sets"]) if "target_sets" in features else []
                covered = sorted(int(x) for panel in panels for x in features["target_sets"][panel][:])
                if panels:
                    assert covered == list(range(len(ids)))
                matrix = csc_matrix((g["data"][:],g["indices"][:],g["indptr"][:]),shape=tuple(g["shape"][:]))
            matrix.sum_duplicates()
            matrix.eliminate_zeros()
            umi = np.asarray(matrix.sum(axis=1)).ravel()
            detection = matrix.getnnz(axis=1)
            out.writerows((name,i,s,int(v),int(n)) for i,s,v,n in zip(ids,symbols,umi,detection))
            metadata.append(dict(sample=name, rows=len(ids), spots=len(barcodes[name]), target_set="|".join(panels)))
    barcode_rows = []
    for m in metadata:
        target = m["sample"]
        if not target.startswith("Targeted"):
            continue
        parent = target.replace("Targeted_", "Parent_", 1).rsplit("_", 1)[0]
        if parent.endswith("_Pan"):
            parent = parent[:-4]
        a, b = barcodes[target], barcodes[parent]
        barcode_rows.append(dict(target=target,parent=parent,target_annotation=m["target_set"],target_spots=len(a),
                                 parent_spots=len(b),shared_barcodes=len(set(a)&set(b)),same_ordered_barcodes=int(a==b)))
    save(HERE / "barcode_pairs.tsv", barcode_rows, barcode_rows[0])


def replay(outdir):
    with gzip.open(HERE / "source_profiles.tsv.gz", "rt", encoding="utf-8") as f:
        by_sample = {}
        for row in csv.DictReader(f,delimiter="\t"):
            by_sample.setdefault(row["sample"], []).append((row["stable_id"],row["symbol"],int(row["umi"]),int(row["detected_spots"])))
    results, discordant = [], []
    for pair in records(HERE/"barcode_pairs.tsv"):
        t, p = pair["target"], pair["parent"]
        assert int(pair["target_spots"]) == int(pair["parent_spots"]) == int(pair["shared_barcodes"])
        for mode in ["deposited_rows","nonzero_before_naming"]:
            axes = []
            lookup = []
            for name in (t,p):
                rows = [r for r in by_sample[name] if mode == "deposited_rows" or r[2] > 0]
                axes.append(Axis([(fid,symbol,None) for fid,symbol,_,_ in rows]))
                lookup.append({r[0]:r[2:] for r in rows})
            stats, events = analyze(*axes)
            uniques = [len(a.by_key) for a in axes]
            results.append(dict(target=t,parent=p,mode=mode,target_rows=uniques[0],parent_rows=uniques[1],
                                shared_barcodes=pair["shared_barcodes"],same_barcode_order=pair["same_ordered_barcodes"],**stats))
            for event in events:
                e = dict(zip(EVENT_FIELDS,event))
                ti, pi = e["left_id"],e["right_id"]
                discordant.append(dict(target=t,mode=mode,event=e["event"],target_key=e["left_key"],parent_key=e["right_key"],
                                       target_id=ti,parent_id=pi,target_umi=lookup[0][ti][0],parent_umi=lookup[1][pi][0],
                                       target_detected_spots=lookup[0][ti][1],parent_detected_spots=lookup[1][pi][1]))
    outdir.mkdir(parents=True,exist_ok=True)
    save(outdir / "pair_results.tsv", results, results[0])
    fields=["target","mode","event","target_key","parent_key","target_id","parent_id","target_umi","parent_umi",
            "target_detected_spots","parent_detected_spots"]
    save(outdir / "events.tsv", discordant, fields)
    print(json.dumps({mode:dict(pairs=len(rows),name_id_discordances=sum(r["different_ids_duplicate"]+r["different_ids_singleton"] for r in rows),
                                shared_ids_different_keys=sum(r["shared_id_different_keys"] for r in rows))
                      for mode in ("deposited_rows","nonzero_before_naming")
                      for rows in [[r for r in results if r["mode"]==mode]]}))


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--extract",type=Path,help="Read 13 original files downloaded by fetch_pairs.py")
    p.add_argument("--output-dir",type=Path,default=HERE)
    a = p.parse_args()
    if a.extract:
        extract(a.extract)
    replay(a.output_dir)
