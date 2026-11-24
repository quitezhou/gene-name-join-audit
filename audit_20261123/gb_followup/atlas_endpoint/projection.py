"""Case-level audit of three author-distributed CellTypist predictors.

python projection.py --rscript "/path/to/Rscript" [--use-extracts]
The large source count matrices and transient query matrices use the OS cache.
"""
import argparse
import csv
import gzip
import hashlib
import io
import json
import pickle
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from scipy.io import mmread
from scipy.special import expit

HERE = Path(__file__).resolve().parent
CASE = HERE.parents[1] / "cohort_case"
sys.path.insert(0, str(CASE))
from cache_entry import cache

COMMIT = "e34ba4698f1b622652b5b0791f1a3e158867e1fa"
MODEL_BLOBS = {"Group": "3809e7739be16f46ded8de522ec8c226129ee6dd",
               "High": "520b8d4bad62027a2e64dbb693327f61db05cecd",
               "Low": "1f3e4b098b75587cc412700f405df69e9366a6fd"}
OBJECTS = ["adipo", "endo", "ly", "mast", "scFAP", "B"]
SOURCE_NAMES = {1: "Massier et al. (#1)", 4: "Massier et al. (#4)/21L005051"}
CANDIDATES = ["AHRR", "ARMCX5-GPRASP2", "ATXN7", "PINX1", "POLR2J3"]
# The overlapping raw-ID intersection was previously verified in these two inputs.
# Keys are the labels used by the three released classifiers.
COMMON = {
    "PINX1.1": ("PINX1", "ENSG00000254093", {1: "PINX1.1", 4: "PINX1"}),
    "ATXN7": ("ATXN7", "ENSG00000163635", {1: "ATXN7.1", 4: "ATXN7"}),
    "POLR2J3": ("POLR2J3", "ENSG00000285437", {1: "POLR2J3.1", 4: "POLR2J3"}),
}


class Record:
    """Inert storage for serialized sklearn attribute dictionaries."""


class ReadModel(pickle.Unpickler):
    def find_class(self, module, name):
        allowed = {("sklearn.linear_model._logistic", "LogisticRegression"): Record,
                   ("sklearn.preprocessing._data", "StandardScaler"): Record,
                   ("numpy.core.multiarray", "_reconstruct"): np.core.multiarray._reconstruct,
                   ("numpy.core.multiarray", "scalar"): np.core.multiarray.scalar,
                   ("numpy", "ndarray"): np.ndarray, ("numpy", "dtype"): np.dtype}
        return allowed[(module, name)]


def read_model(name):
    folder = HERE / "models"
    folder.mkdir(exist_ok=True)
    path = folder / f"Cluster_{name}.pkl"
    url = f"https://raw.githubusercontent.com/lmassier/hWAT_singlecell/{COMMIT}/CellTypist/{path.name}"
    if not path.exists():
        r = requests.get(url, timeout=45)
        r.raise_for_status()
        path.write_bytes(r.content)
    data = path.read_bytes()
    blob = hashlib.sha1(f"blob {len(data)}\x00".encode() + data).hexdigest()
    assert blob == MODEL_BLOBS[name]
    model = ReadModel(io.BytesIO(data)).load()
    lr, sc = model["Model"], model["Scaler_"]
    assert len(set(lr.features)) == len(lr.features) == lr.coef_.shape[1] == len(sc.scale_)
    assert np.all(sc.scale_ > 0) and len(lr.classes_) == lr.coef_.shape[0]
    return lr, sc, dict(name=name, url=url, sha256=hashlib.sha256(data).hexdigest(), git_blob=blob,
                        date=model["description"]["date"], trained_sklearn=lr._sklearn_version,
                        classes=len(lr.classes_), features=len(lr.features))


def prepare_sources(inv, temp):
    bars = {}
    frequency = Counter()
    for name, record in inv.items():
        if name.startswith("Massier") and name.endswith("/barcodes.tsv.gz"):
            with gzip.open(cache(record), "rt") as f:
                b = [s.strip() for s in f]
            bars[name.rsplit("/", 1)[0]] = b
            frequency.update(b)
    rows, source_info = [], {}
    for k, name in SOURCE_NAMES.items():
        kept = {b: i for i, b in enumerate(bars[name]) if frequency[b] == 1}
        rows.extend({"barcode": b, "source_number": k} for b in kept)
        ft = inv[name + "/features.tsv.gz"]
        fpath = CASE / "source" / (ft["file_id"] + ".tsv.gz")
        assert hashlib.sha256(fpath.read_bytes()).hexdigest() == ft["sha256"]
        feats = pd.read_csv(fpath, sep="\t", header=None)
        with gzip.open(cache(inv[name + "/matrix.mtx.gz"]), "rb") as f:
            m = mmread(f).tocsr()
        assert m.shape == (len(feats), len(bars[name]))
        vec = {}
        for _, (symbol, stable, _) in COMMON.items():
            group = feats[feats[1] == symbol]
            assert list(group[0]).count(stable) == 1
            vec[stable] = m[group.index[group[0] == stable][0]].toarray()[0]
        source_info[k] = (kept, vec)
    pd.DataFrame(rows).to_csv(temp / "barcodes.tsv", index=False, sep="\t")
    return source_info


def top(scores):
    order = np.argsort(-scores, axis=1, kind="stable")
    i = np.arange(scores.shape[0])
    return order[:, 0], scores[i, order[:, 0]] - scores[i, order[:, 1]]


def run():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rscript", default="Rscript")
    parser.add_argument("--use-extracts", action="store_true")
    args = parser.parse_args()
    with (CASE / "deposit_inventory.tsv").open(encoding="utf8") as f:
        inv = {r["path"]: r for r in csv.DictReader(f, delimiter="\t")}
    temp = Path(tempfile.gettempdir()) / "massier_endpoint_projection"
    temp.mkdir(exist_ok=True)
    model = {k: read_model(k) for k in MODEL_BLOBS}
    targets = {g + suffix for g in CANDIDATES for suffix in ["", ".1", "-1"]}
    keyset = set.union(*(set(lr.features) for lr, _, _ in model.values()))
    assert keyset & targets == set(COMMON)
    (temp / "keys.txt").write_text("\n".join(sorted(keyset | targets)) + "\n")
    raw = prepare_sources(inv, temp)
    summary, per_cell, audits, weights = [], [], [], []
    seen_barcodes = set()
    for m, (lr, sc, info) in model.items():
        for j in np.flatnonzero(np.isin(lr.features, list(COMMON))):
            for i, celltype in enumerate(lr.classes_):
                weights.append(dict(model=m, feature=lr.features[j], celltype=celltype,
                                    coefficient=lr.coef_[i, j], scale=sc.scale_[j],
                                    mean=sc.mean_[j], with_mean=sc.with_mean))
    for obj in OBJECTS:
        stem = temp / obj
        if not args.use_extracts:
            rds = next(v for k, v in inv.items() if k.endswith("/" + obj + ".rds"))
            subprocess.run([args.rscript, str(HERE / "extract_queries.R"), str(cache(rds)),
                            str(temp / "keys.txt"), str(temp / "barcodes.tsv"), str(stem)], check=True)
        keys = Path(f"{stem}_keys.txt").read_text().splitlines()
        ix = {s: j for j, s in enumerate(keys)}
        counts = mmread(f"{stem}.mtx").tocsr()
        obs = pd.read_csv(f"{stem}_obs.tsv", sep="\t")
        assert counts.shape == (len(keys), len(obs)) and obs.barcode.is_unique and min(obs.total)>0
        assert seen_barcodes.isdisjoint(obs.barcode)
        seen_barcodes.update(obs.barcode)
        for k, subset in obs.groupby("source"):
            cols = subset.index.to_numpy()
            source_idx, source_vectors = raw[k]
            cols_raw = [source_idx[b] for b in subset.barcode]
            recovered = {}
            for key, (label, gid, q) in COMMON.items():
                anchor = q[k]
                values = counts[ix[anchor], cols].toarray().ravel() if anchor in ix else np.zeros(len(cols))
                support = source_vectors[gid][cols_raw] > 0
                agree = np.array_equal(values > 0, support)
                audits.append(dict(object=obj, source=k, model_key=key, stable_id=gid, assay_row=anchor,
                                   row_available=anchor in ix, cells=len(cols), raw_detected=int(support.sum()),
                                   deposited_detected=int((values > 0).sum()), equal_support=agree))
                recovered[key] = values if agree else None
            for m, (lr, sc, info) in model.items():
                mapping = {s: j for j, s in enumerate(lr.features)}
                names = [g for g in keys if g in mapping]
                model_idx = [mapping[g] for g in names]
                x = counts[[ix[g] for g in names]][:, cols].T.toarray()
                x = np.log1p(x * (1e4 / subset.total.to_numpy())[:, None])
                mu = sc.mean_[model_idx] if sc.with_mean else np.zeros(len(names))
                std = sc.scale_[model_idx]
                z = np.minimum((x - mu) / std, 10.)
                coef = lr.coef_[:, model_idx]
                scores = z @ coef.T + lr.intercept_
                target_keys = set(COMMON) & set(lr.features)
                if not target_keys <= set(names) or any(recovered[key] is None for key in target_keys):
                    summary.append(dict(object=obj, source=k, model=m, cells=len(cols),
                                        status="nonadjudicable_missing_row_or_changed_detection"))
                    continue
                fixed = z.copy()
                for key in target_keys:
                    j = names.index(key)
                    normalized = np.log1p(recovered[key] * 1e4 / subset.total)
                    fixed[:,j] = np.minimum((normalized - mu[j]) / std[j], 10.)
                    assert np.all(np.isfinite(fixed[:,j]))
                scores_fixed = fixed @ coef.T + lr.intercept_
                np.savez_compressed(temp / f"{obj}_{k}_{m}_scores.npz",
                                    baseline=scores, aligned=scores_fixed)
                ti = [names.index(g) for g in sorted(target_keys)]
                removed = scores - z[:,ti] @ coef[:,ti].T
                whole_delta = scores_fixed - scores
                assert np.allclose(whole_delta, (fixed[:,ti] - z[:,ti]) @ coef[:,ti].T, atol=1e-9)
                orig, margin = top(scores)
                corr, margin_c = top(scores_fixed)
                omit, margin_o = top(removed)
                affected = np.any(np.abs(fixed[:,ti] - z[:,ti]) > 1e-9, axis=1)
                summary.append(dict(object=obj, source=k, model=m, cells=len(cols), status="adjudicated",
                                    features_overlapping=len(names), model_features=len(lr.features),
                                    addressed="|".join(sorted(target_keys)), measurement_changed=int(affected.sum()),
                                    top_changed=int((orig!=corr).sum()), top_changed_omit=int((orig!=omit).sum()),
                                    max_abs_class_score_shift=np.abs(whole_delta).max(),
                                    max_abs_margin_shift=np.abs(margin_c-margin).max(),
                                    max_abs_class_probability_shift=np.abs(expit(scores_fixed)-expit(scores)).max()))
                for j in np.flatnonzero(affected | (orig!=omit)):
                    per_cell.append(dict(object=obj, source=k, barcode=subset.iloc[j].barcode, model=m,
                                         baseline_label=lr.classes_[orig[j]], shared_ID_label=lr.classes_[corr[j]],
                                         omitted_label=lr.classes_[omit[j]], baseline_margin=margin[j],
                                         shared_ID_margin=margin_c[j], omitted_margin=margin_o[j], changed_input=bool(affected[j]),
                                         max_abs_logit_shift=np.max(np.abs(whole_delta[j]))))
    for name, vals in [("model_weights", weights), ("query_identity_checks", audits),
                       ("score_sensitivity", summary), ("cell_sensitivity", per_cell)]:
        pd.DataFrame(vals).to_csv(HERE / (name + ".tsv"), sep="\t", index=False)
    scored = pd.DataFrame(summary).query("status == 'adjudicated'")
    aggregate = scored.groupby("model")[["cells", "measurement_changed", "top_changed", "top_changed_omit"]].sum()
    aggregate.to_csv(HERE / "paired_summary.tsv", sep="\t")
    (HERE / "models.json").write_text(json.dumps([item[2] for item in model.values()], indent=2) + "\n")


if __name__ == "__main__":
    run()
