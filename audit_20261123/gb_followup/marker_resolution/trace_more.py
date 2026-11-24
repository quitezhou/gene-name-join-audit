"""Trace source cells using the target-blind 60-gene panel and the full raw gene axis."""
from collections import Counter, defaultdict
import csv
import gzip
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.io import mmread

HERE = Path(__file__).resolve().parent
PREVIOUS = HERE.parent / "hwat_model"
sys.path.insert(0, str(PREVIOUS))
from tcell_published_marker import BASE, ID, PUBLISHED, STUDY, cache, feature_records, name_table, support_match


def main():
    cells = pd.read_csv(PREVIOUS / "source/ly_rows.tsv.gz", sep="\t")
    ctrls = pd.read_csv(HERE / "control_detection.tsv.gz", sep="\t").set_index("cell").loc[cells.cell]
    parts = cells.cell.str.split("_", n=1, expand=True)
    classes = cells.source + "|" + parts[1]
    strata = cells.source + "|" + cells.cell.str.extract(r"((?:_\d+)+)$")[0].fillna("unknown")
    old = pd.read_csv(PREVIOUS / "tcell_published_marker_id_trace.tsv", sep="\t").set_index("sample")
    inv = {r["path"]: r for r in csv.DictReader((STUDY / "deposit_inventory.tsv").open(), delimiter="\t")}
    panels = {k: rows for k, (_, rows) in feature_records().items() if k.startswith("Massier")}
    barcodes = {}
    for sample in panels:
        with gzip.open(cache(inv[sample + "/barcodes.tsv.gz"]), "rt") as stream:
            bs = [s.strip() for s in stream]
        assert len(bs) == len(set(bs))
        barcodes[sample] = {b: i for i, b in enumerate(bs)}
    frequency = Counter(b for panel in barcodes.values() for b in panel)
    for sample, rows in panels.items():
        source_idx = barcodes[sample]
        mapped = parts[0].map(source_idx)
        original_class = old.at[sample, "source_suffix"]
        eligible = classes.eq(original_class) & mapped.notna() & ~cells.cell.duplicated(keep=False)
        candidates = np.flatnonzero(eligible)
        ids = [row[0] for row in rows]
        leaf, freq = name_table(rows)
        focus = [g for g in ctrls if g in leaf and freq[g] == 1]
        assert len(focus) >= 40
        with gzip.open(cache(inv[sample + "/matrix.mtx.gz"]), "rb") as stream:
            arr = mmread(stream).tocsc()
        ix = mapped.iloc[candidates].astype(int).to_numpy()
        raw = arr[:, ix].tocsr()
        raw.eliminate_zeros()
        raw.sort_indices()
        ctr_rows = [ids.index(leaf[g][0]) for g in focus]
        detections = raw[ctr_rows].toarray().T > 0
        obj_control = ctrls.iloc[candidates][focus].to_numpy() > 0
        unique_bar = parts[0].iloc[candidates].map(frequency).eq(1).to_numpy()
        matches_control = np.all(detections == obj_control, axis=1)
        yield sample, raw, candidates, ids, leaf, unique_bar, detections, obj_control, matches_control, ctr_rows, focus, cells, strata, old


def run():
    info = []
    per_cell = defaultdict(list)
    for sample, raw, ix, ids, name_mapping, unique, det, obj, matches, row_ids, focus, cells, strata, old in main():
        key, = [key for key, (gid, _) in name_mapping.items() if gid == ID]
        inferred = raw[ids.index(ID)].toarray().ravel() > 0
        observed = cells["counts_" + key].iloc[ix].gt(0).to_numpy()
        focus_matches = []
        # Every control gene with three or more detections must be supported
        # at the expected positions; only unique raw-row fingerprints count.
        for j, g in enumerate(focus):
            pos = obj[:, j]
            if pos.sum() < 3:
                continue
            hit = support_match(raw, pos)
            if len(hit) == 1 and hit[0] == row_ids[j]:
                focus_matches.append(g)
        target_hit = support_match(raw, observed)
        target_pass = observed.sum() >= 3 and len(target_hit) == 1 and ids[target_hit[0]] == ID
        negative = np.all(np.roll(det,1,axis=0) == obj,axis=1).sum()
        broad = bool(matches.all() and len(focus_matches) >= 2 and target_pass
                     and negative < .1 * len(ix) and len(ix) >= 20)
        exact = unique & matches & (det.sum(axis=1)>=3)
        control_ok = []
        target_ok = False
        shifted_matches = None
        if exact.sum() >= 20:
            subset = raw[:,exact].tocsr()
            for j,g in enumerate(focus):
                pos = obj[exact,j]
                if pos.sum() < 3:
                    continue
                hits = support_match(subset, pos)
                if len(hits) == 1 and hits[0] == row_ids[j]:
                    control_ok.append(g)
            target_hits = support_match(subset, observed[exact])
            target_ok = bool(observed[exact].sum() >= 3 and len(target_hits) == 1 and ids[target_hits[0]] == ID)
            shifted_matches = int(np.all(np.roll(det[exact],1,axis=0)==obj[exact],axis=1).sum())
        cell_arm = bool(len(control_ok) >= 2 and target_ok and shifted_matches < .1 * exact.sum())
        verified_before = bool(old.at[sample,"include"])
        assert int(unique.sum()) == old.at[sample,"barcode_pairs"]
        for j, cell in enumerate(ix):
            before = verified_before and unique[j]
            if not before and not ((broad or cell_arm) and exact[j]):
                continue
            per_cell[int(cell)].append((sample, bool(inferred[j]), bool(before), bool(unique[j])))
        info.append(dict(sample=sample, source_suffix=old.at[sample,"source_suffix"],
                         cells=len(ix), unique_barcodes=int(unique.sum()), control_genes=len(focus),
                         control_positive_genes_with_unique_pattern=len(focus_matches),
                         cells_with_exact_control_vector=int(matches.sum()), shifted_barcode_neg_match=int(negative),
                         target_detections=int(observed.sum()), target_pattern_unique_and_correct=target_pass,
                         complete_stratum_verified=broad, original_accepted=verified_before,
                         exact_subset_cells=int(exact.sum()), exact_subset_controls=len(control_ok),
                         exact_subset_shifted_matches=shifted_matches, exact_subset_target_verified=target_ok,
                         cell_subset_accepted=cell_arm))
    baseline = cells["counts_" + PUBLISHED].gt(0).to_numpy()
    out = pd.DataFrame(dict(cell=cells.cell, source_stratum=strata, in_cluster=cells.cluster.eq(3).astype(int),
                            bare=cells["counts_"+BASE].gt(0).astype(int), printed=baseline.astype(int),
                            before=0, resolved=-1, source_library=""))
    for cell, options in per_cell.items():
        assert len(options) == 1 or all(not v[2] for v in options)
        if len(options) != 1:
            continue
        sample, val, before, unique = options[0]
        # The barcode was previously classed ambiguous if present in another
        # raw library. A compatible *other* source could fail our filters.
        # Exclude these outright even if only one source passed every gate.
        if not unique:
            continue
        out.loc[cell, ["resolved","before","source_library"]] = [int(val), int(before), sample]
    assert out.before.sum() == 1375
    pd.DataFrame(info).to_csv(HERE/"expanded_library_checks.tsv", index=False, sep="\t")
    out.to_csv(HERE/"cell_detection_trace.tsv.gz", index=False, sep="\t",
               compression={"method":"gzip","mtime":0})


if __name__ == "__main__":
    run()
