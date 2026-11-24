"""Trace the published T-cell marker's stable ID on source-bound barcodes.

An orig.ident string can differ from the deposited cohort number. Require a
dominant object barcode suffix per library and globally unique raw barcodes,
then verify detection patterns against every original feature and unrelated
raw/object controls. The conservative repair edits only qualified cells.
"""
import csv
import gzip
from collections import Counter
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.io import mmread

from network_edges import STUDY, name_table
from cache_entry import cache
from derive_features import feature_records

HERE = Path(__file__).resolve().parent
ID = "ENSG00000285437"
BASE = "POLR2J3"
PUBLISHED = "POLR2J3.1"


def support_match(matrix, target):
    active = np.flatnonzero(target)
    degrees = np.diff(matrix.indptr)
    overlap = matrix[:, active].getnnz(axis=1)
    matches = np.flatnonzero((overlap == len(active)) & (degrees == len(active)))
    return matches


def main():
    inv = {r["path"]: r for r in csv.DictReader(
        (STUDY / "deposit_inventory.tsv").open(), delimiter="\t")}
    panels = feature_records()
    cells = pd.read_csv(HERE / "source/ly_rows.tsv.gz", sep="\t")
    token = cells.cell.str.split("_", n=1, expand=True)
    bars = {}
    for sample in panels:
        if sample.startswith("Massier"):
            with gzip.open(cache(inv[sample + "/barcodes.tsv.gz"]), "rt") as src:
                bars[sample] = [line.strip() for line in src]
    freq = Counter(b for bs in bars.values() for b in bs)
    original = cells["counts_" + PUBLISHED].gt(0).to_numpy()
    repaired = original.copy()
    used = np.zeros(len(cells), dtype=bool)
    cluster = cells.cluster.eq(3).to_numpy()
    traces = []
    for sample, bs in bars.items():
        index = {b: i for i, b in enumerate(bs) if freq[b] == 1}
        mapping = token[0].map(index)
        # Choose a metadata/suffix class using only barcode overlap. Discard
        # chance matches in other source classes and repeated object barcodes.
        cls = cells.source + "|" + token[1]
        counts = cls[mapping.notna()].value_counts()
        if counts.empty:
            continue
        chosen = counts.index[0]
        candidates = np.flatnonzero((mapping.notna() & cls.eq(chosen)).to_numpy())
        candidates = candidates[~token.iloc[candidates, 0].duplicated(keep=False).to_numpy()]
        rows = panels[sample][1]
        names, _ = name_table(rows)
        row_keys = [k for k, (gid, _) in names.items() if gid == ID]
        assert len(row_keys) == 1
        key = row_keys[0]
        assert key in (PUBLISHED, BASE)
        with gzip.open(cache(inv[sample + "/matrix.mtx.gz"]), "rb") as src:
            raw = mmread(src).tocsc()
        indices = mapping.iloc[candidates].astype(int).to_numpy()
        selected = raw[:, indices].tocsr()
        selected.eliminate_zeros()
        selected.sort_indices()
        feature_ids = [r[0] for r in rows]
        target = cells["counts_" + key].iloc[candidates].to_numpy()
        matches = support_match(selected, target > 0)
        id_row = feature_ids.index(ID)
        raw_target = selected[id_row].toarray().ravel()
        unique_id = len(matches) == 1 and matches[0] == id_row
        # Two unrelated identity controls are sufficient only if both have
        # informative, globally unique support and exactly follow the source
        # feature table. Keep the actual counts and all failed arms.
        controls = []
        for control, (cid, family) in names.items():
            if family not in {"ATXN7", "PINX1", "ARMCX5-GPRASP2", "AHRR"}:
                continue
            field = "counts_" + control
            if field not in cells:
                continue
            det = cells[field].iloc[candidates].to_numpy() > 0
            if det.sum() < 3:
                continue
            hits = support_match(selected, det)
            if len(hits) == 1 and feature_ids[hits[0]] == cid:
                controls.append(control)
        passes = (len(candidates) >= 20 and len(candidates) / counts.sum() >= .8
                  and unique_id and (target > 0).sum() >= 3 and len(controls) >= 2)
        if passes:
            assert not used[candidates].any()
            used[candidates] = True
            repaired[candidates] = raw_target > 0
        traces.append(dict(sample=sample, source_suffix=chosen, barcode_pairs=len(candidates),
                           all_classes_barcode_pairs=int(counts.sum()), raw_rows=len(rows),
                           native_key=key, informative_object_cells=int((target > 0).sum()),
                           support_matches=len(matches), unique_expected_id=unique_id,
                           raw_values_exact=np.array_equal(raw_target, target),
                           control_support_matches="|".join(controls), include=passes,
                           target_cluster_cells=int(cluster[candidates].sum()),
                           printed_key_positive_in=int(original[candidates][cluster[candidates]].sum()),
                           printed_key_positive_out=int(original[candidates][~cluster[candidates]].sum()),
                           id_positive_in=int((raw_target[cluster[candidates]] > 0).sum()),
                           id_positive_out=int((raw_target[~cluster[candidates]] > 0).sum())))
        print(sample, len(candidates), key, unique_id, controls, passes, flush=True)
    pd.DataFrame(traces).to_csv(HERE / "tcell_published_marker_id_trace.tsv", sep="\t", index=False)
    output = []
    for label, v, selection in [
        ("released_object_printed_key", original, np.ones(len(cells), dtype=bool)),
        ("released_object_source_verified_partial_id_repair", repaired, np.ones(len(cells), dtype=bool)),
        ("verified_cells_printed_key", original, used),
        ("verified_cells_id_repair", repaired, used),
        ("observed_dominant_source_printed_key", original, cells.source.eq("Massier et al., #3").to_numpy()),
    ]:
        a, b = selection & cluster, selection & ~cluster
        output.append(dict(arm=label, in_n=int(a.sum()), out_n=int(b.sum()),
                           in_positive=int(v[a].sum()), out_positive=int(v[b].sum()),
                           in_fraction=v[a].mean(), out_fraction=v[b].mean(),
                           difference_pp=100*(v[a].mean() - v[b].mean())))
    print(pd.DataFrame(output).to_string(index=False), flush=True)
    pd.DataFrame(output).to_csv(HERE / "tcell_published_marker_detection.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
