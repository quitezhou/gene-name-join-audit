"""Independently recompute QC pseudobulks from one GEO MEX per observable design stratum."""
import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
from collections import Counter

import numpy as np
import pandas as pd
from scipy.io import mmread

HERE = Path(__file__).resolve().parent
CASE = HERE.parents[1] / "sentinel_case"
PROJECT = CASE.parent.parent


def main(cache):
    results = pd.read_csv(HERE / "all_results.tsv.gz", sep="\t", index_col="gene")
    renamed = results[results.match == "one_release_label_same_id"]
    archive = np.load(CASE / "pseudobulk_counts.npz")
    manifest = json.loads((CASE / "matrix_manifest.json").read_text())
    records = []
    summaries = {}
    for sample in ["TRNA_1", "NEADT_1", "NEADT_12"]:
        rec = next(r for r in manifest if r["sample"] == sample)
        acc = rec["accession"]
        name = f"{acc}_{sample}_matrix.mtx.gz"
        path = cache / name
        if not path.is_file():
            path = PROJECT / "panel_sources" / name
        assert hashlib.sha256(path.read_bytes()).hexdigest() == rec["sha256"]
        with gzip.open(path, "rb") as stream:
            matrix = mmread(stream).tocsc()
        detected_spots = np.asarray((matrix > 0).sum(axis=1)).ravel()
        matrix = matrix.multiply((detected_spots >= 5)[:, None]).tocsc()
        matrix.eliminate_zeros()
        keep = np.asarray(matrix.sum(axis=0)).ravel() >= 500
        counts = np.asarray(matrix[:, keep].sum(axis=1)).ravel()
        assert np.array_equal(counts, archive[acc + "_qc"])
        features_path = PROJECT / "panel_sources/all_features" / f"{acc}_{sample}_features.tsv.gz"
        with gzip.open(features_path, "rt") as stream:
            features = list(csv.reader(stream, delimiter="\t"))
        names = Counter(r[1] for r in features)
        by_id = {r[0]: i for i, r in enumerate(features)}
        release = "new" if len(by_id) == 36601 else "old"
        for label, row in renamed.iterrows():
            stable_idx = by_id[row.anchor]
            raw_name = features[stable_idx][1]
            assert names[raw_name] == 1 and raw_name == row[f"anchor_{release}_label"]
            assert raw_name != row["anchor_old_label" if release == "new" else "anchor_new_label"]
            if label != raw_name:
                assert names[label] == 0
            records.append(dict(sample=sample, gene=label, id=row.anchor, raw_name=raw_name,
                                local=int(counts[stable_idx]) if raw_name == label else 0,
                                shared_id=int(counts[stable_idx])))
        summaries[sample] = dict(sha256=rec["sha256"], all_gene_counts_equal=True,
                                 feature_rows=len(counts), retained_spots=int(keep.sum()))
    pd.DataFrame(records).to_csv(HERE / "independent_MEX_trace.tsv", sep="\t", index=False)
    (HERE / "MEX_verification.json").write_text(json.dumps(summaries, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cache", type=Path)
    main(parser.parse_args().cache)
