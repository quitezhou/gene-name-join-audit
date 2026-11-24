"""Counterfactual Tangram naming audit on deposited skin count matrices.

The original study is NOT asserted to have used this workflow. Per-gene sums
preserve exactly min_cells=1 survival for nonnegative counts, not other
filtering, density priors, spatial neighbors, or model optimization.
"""
from __future__ import annotations
import argparse
import ast
import gzip
import hashlib
import json
import logging
from pathlib import Path
import warnings

import anndata
from anndata.utils import make_index_unique
import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

BASE = Path(__file__).resolve().parent
MODES = ("filter_then_unique", "unique_before_filter", "stable_id")


def fingerprint(rows):
    return hashlib.sha256("\n".join(f"{r['stable_id']}\t{r['symbol']}\t{r['feature_type']}" for r in rows).encode()).hexdigest()


def spatial_summary(path):
    with h5py.File(path, "r") as h:
        m = h["matrix"]
        n_genes, n_spots = map(int, m["shape"][:])
        data, indices = m["data"][:], m["indices"][:]
        assert np.issubdtype(data.dtype, np.integer) and np.all(data >= 0)
        total = np.zeros(n_genes, dtype=np.int64)
        positive = np.zeros(n_genes, dtype=np.int64)
        np.add.at(total, indices, data)
        np.add.at(positive, indices, data > 0)
        cols = [m["features"][key].asstr()[:].tolist() for key in ("id", "name", "feature_type")]
    rows = [dict(stable_id=g, symbol=s, feature_type=t, total_count=int(c), positive_entries=int(p))
            for g,s,t,c,p in zip(*cols, total, positive)]
    assert len(rows) == n_genes and len({r["stable_id"] for r in rows}) == n_genes
    return dict(sample=path.name.split("_")[0], series="GSE241124", n_genes=n_genes,
                n_cells=n_spots, total_counts=int(total.sum()), nnz=len(data),
                file=str(path.relative_to(BASE)), file_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                feature_rows=rows)


def maps(sample, mode):
    rows = sample["feature_rows"]
    symbols = pd.Index([r["symbol"] for r in rows])
    if mode == "unique_before_filter":
        names = make_index_unique(symbols)
    elif mode == "stable_id":
        names = pd.Index([r["stable_id"] for r in rows])
    else:
        names = symbols
    keep = np.array([r["total_count"] > 0 for r in rows])
    # pp_adatas genes=None takes the reference candidates after filtering but
    # before lowercasing/unique repair. It does NOT use all repaired suffixes.
    candidates = pd.Index([str(s).lower() for s in names[keep]])
    final = make_index_unique(candidates)
    records = {key:r for key,r in zip(final, (r for r,k in zip(rows,keep) if k))}
    return records, set(candidates)


def compare(ref, sp, mode, cache):
    rm, candidates = cache[(ref["sample"], mode)]
    sm, _ = cache[(sp["sample"], mode)]
    overlap = sorted(rm.keys() & sm.keys())
    training = sorted(set(overlap) & candidates)
    conflicts = []
    for key in overlap:
        a,b = rm[key], sm[key]
        if a["stable_id"] != b["stable_id"]:
            conflicts.append(dict(key=key, reference_id=a["stable_id"], spatial_id=b["stable_id"],
                                  reference_symbol=a["symbol"], spatial_symbol=b["symbol"],
                                  reference_count=a["total_count"], spatial_count=b["total_count"],
                                  enters_default_training=key in candidates))
    train_conflicts = [r for r in conflicts if r["enters_default_training"]]
    rid = {r["stable_id"]: k for k,r in rm.items()}
    sid = {r["stable_id"]: k for k,r in sm.items()}
    split = [g for g in rid.keys() & sid.keys() if rid[g] != sid[g]]
    return dict(reference=ref["sample"], spatial=sp["sample"], mode=mode,
                training_genes=len(training), overlap_genes=len(overlap),
                overlap_stable_id_conflicts=len(conflicts), training_stable_id_conflicts=len(train_conflicts),
                same_id_different_key=len(split), split_ids=sorted(split), conflicts=conflicts,
                conflict_reference_counts=sum(x["reference_count"] for x in train_conflicts),
                conflict_spatial_counts=sum(x["spatial_count"] for x in train_conflicts),
                reference_total_counts=ref["total_counts"], spatial_total_counts=sp["total_counts"])


def validate_source(ref, sp, mode, expected, cache):
    from probe_tangram_pp_adatas import source_only_scanpy, COMMIT, URL
    sc, hashes = source_only_scanpy()
    raw = (BASE / "tangram_mapping_utils_pinned.py").read_bytes()
    fn = next(n for n in ast.parse(raw.decode()).body if isinstance(n, ast.FunctionDef) and n.name == "pp_adatas")
    scope = {"np": np, "pd": pd, "sc": sc, "logging": logging}
    exec(compile(ast.Module(body=[fn], type_ignores=[]), URL, "exec"), scope)
    def obj(sample):
        rows = sample["feature_rows"]
        names = pd.Index([r["stable_id"] if mode == "stable_id" else r["symbol"] for r in rows])
        if mode == "unique_before_filter":
            names = make_index_unique(names)
        var = pd.DataFrame({"stable_id": [r["stable_id"] for r in rows]}, index=names)
        return anndata.AnnData(csr_matrix([[r["total_count"] for r in rows]], dtype=np.float64), var=var)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Variable names are not unique")
        a,b = obj(ref), obj(sp)
        scope["pp_adatas"](a,b)
    rm, candidates = cache[(ref["sample"],mode)]
    sm,_ = cache[(sp["sample"],mode)]
    assert list(a.var_names) == list(rm) and list(b.var_names) == list(sm)
    assert set(a.uns["training_genes"]) == candidates & rm.keys() & sm.keys()
    assert set(a.uns["overlap_genes"]) == rm.keys() & sm.keys()
    wrong = sum(a.var.loc[k,"stable_id"] != b.var.loc[k,"stable_id"] for k in a.uns["training_genes"])
    assert wrong == expected["training_stable_id_conflicts"]
    return dict(reference=ref["sample"], spatial=sp["sample"], mode=mode,
                training_conflicts=int(wrong), naming_filter_and_training_set_assertions="passed",
                tangram_commit=COMMIT, tangram_source_sha256=hashlib.sha256(raw).hexdigest(),
                scanpy_source_sha256=hashes)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reference-dir", type=Path, default=BASE/"tangram_real_skin/reference")
    ap.add_argument("--allow-partial", action="store_true")
    args = ap.parse_args()
    refs = []
    for path in sorted(args.reference_dir.glob("*_count_summary.json.gz")):
        with gzip.open(path,"rt",encoding="utf-8") as f:
            refs.append(json.load(f))
    if not args.allow_partial:
        assert len(refs) == 12, len(refs)
    assert refs
    spatial = [spatial_summary(p) for p in sorted((BASE/"cohort_direct_features").glob("*.h5"))
               if p.name.startswith(("GSM77170", "GSM82384"))]
    assert len(spatial) == 16
    hashes = {fingerprint(x["feature_rows"]) for x in refs+spatial}
    assert len(hashes) == 1, hashes
    pooled_rows = [dict(r) for r in refs[0]["feature_rows"]]
    for i,r in enumerate(pooled_rows):
        r["total_count"] = sum(x["feature_rows"][i]["total_count"] for x in refs)
        r["positive_entries"] = sum(x["feature_rows"][i]["positive_entries"] for x in refs)
    pooled = dict(sample=f"POOLED_{len(refs)}_REFERENCES", series="GSE241132", n_genes=len(pooled_rows),
                  n_cells=sum(r["n_cells"] for r in refs), total_counts=sum(r["total_counts"] for r in refs),
                  feature_rows=pooled_rows)
    allrefs = refs+[pooled]
    cache = {(x["sample"],m):maps(x,m) for x in allrefs+spatial for m in MODES}
    pairs = [compare(r,s,m,cache) for m in MODES for r in allrefs for s in spatial]
    source_tests = []
    # Prespecified maximum-conflict individual and pooled pairs, each in all
    # three modes. Also check one individual zero-conflict pair if available.
    late = [p for p in pairs if p["mode"] == MODES[0]]
    selected = []
    for is_pool in (False,True):
        sub = [p for p in late if p["reference"].startswith("POOLED") == is_pool]
        selected.append(max(sub,key=lambda p:p["training_stable_id_conflicts"]))
    zero = next((p for p in late if p["training_stable_id_conflicts"] == 0 and not p["reference"].startswith("POOLED")),None)
    if zero and (zero["reference"],zero["spatial"]) not in {(p["reference"],p["spatial"]) for p in selected}:
        selected.append(zero)
    byname = {x["sample"]:x for x in allrefs+spatial}
    for p in selected:
        for m in MODES:
            expected = next(x for x in pairs if x["reference"] == p["reference"] and x["spatial"] == p["spatial"] and x["mode"] == m)
            source_tests.append(validate_source(byname[p["reference"]],byname[p["spatial"]],m,expected,cache))
    summary = []
    for mode in MODES:
        for pooled_flag in (False,True):
            group = [p for p in pairs if p["mode"] == mode and p["reference"].startswith("POOLED") == pooled_flag]
            summary.append(dict(mode=mode,reference_scope="pooled" if pooled_flag else "individual",
                                sample_pairs=len(group), positive_training_pairs=sum(p["training_stable_id_conflicts"]>0 for p in group),
                                max_training_conflicts=max(p["training_stable_id_conflicts"] for p in group),
                                training_conflict_keys=sorted({c["key"] for p in group for c in p["conflicts"] if c["enters_default_training"]}),
                                max_overlap_conflicts=max(p["overlap_stable_id_conflicts"] for p in group)))
    result = dict(date="2026-09-19",scope="Counterfactual workflow experiment on real deposited counts; NOT adjudication of the original paper's actual processing; no optimizer or biological endpoint run.",
                  filter_equivalence="All counts are nonnegative. Sum>0 is exactly equivalent to surviving min_cells=1. One-row aggregation is used ONLY to validate naming/filtering/feature selection, not density or optimization.",
                  interpretation="Different stable IDs are identity-conflict candidates; biological equivalence/locus checks remain required. Pair counts are dependent and are not study prevalence. genes=None, default lowercasing, AnnData default - suffix. Naming modes also change the genes=None candidate universe; these are identity controls, not a fixed-feature downstream effect experiment.",
                  complete_reference_panel=len(refs)==12, raw_ordered_feature_sha256=next(iter(hashes)),
                  n_reference_samples=len(refs), n_spatial_samples=len(spatial),
                  sources=[{k:v for k,v in x.items() if k != "feature_rows"} for x in refs+spatial],
                  summary=summary, source_function_checks=source_tests,pairs=pairs)
    path=BASE/"tangram_skin_real_counts_result.json"
    path.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({"output":str(path),"complete":result["complete_reference_panel"],"summary":summary,"source_checks":len(source_tests)},ensure_ascii=False,indent=2))


if __name__ == "__main__":
    main()
