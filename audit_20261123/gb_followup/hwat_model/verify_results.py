"""Check cross-file denominators and manuscript-facing numbers."""
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

from network_edges import cytab

ROOT = Path(__file__).resolve().parent


def table(name):
    return pd.read_csv(ROOT / name, sep="\t")


def main():
    marker = table("MOESM5_hits.tsv")
    assert len(marker) == 39 and (marker.p_adj < 0.05).sum() == 29
    published = marker[(marker.sheet.str.strip() == "T cells") & (marker.row == 6761)].iloc[0]
    assert published.gene == "POLR2J3.1" and published.cluster == 3
    det = table("marker_detection_validation.tsv")
    assert det.same_cluster_count.sum() == 15
    assert det.both_percentages_match.eq(True).sum() == 1
    tr = table("tcell_published_marker_id_trace.tsv")
    keep = tr[tr.include]
    assert len(keep) == 4 and keep.barcode_pairs.sum() == 1375
    assert keep.unique_expected_id.all()
    baseline, fixed = table("tcell_published_marker_detection.tsv").iloc[:2].to_dict("records")
    assert (baseline["in_n"], baseline["out_n"]) == (1154, 9240)
    assert (baseline["in_positive"], baseline["out_positive"]) == (151, 233)
    assert (round(baseline["in_fraction"], 3), round(baseline["out_fraction"], 3)) == (published.pct1, published.pct2)
    for group in ("in", "out"):
        delta = keep[f"id_positive_{group}"].sum() - keep[f"printed_key_positive_{group}"].sum()
        assert baseline[f"{group}_positive"] + delta == fixed[f"{group}_positive"]
    assert (fixed["in_positive"], fixed["out_positive"]) == (163, 487)
    edges = table("fig1b_tested_edges.tsv")
    with zipfile.ZipFile(ROOT / "source/Figure_1b.cys") as archive:
        assert len(cytab(archive, "Edge")) == 1303
    assert len(edges) == 277 and edges.weight_and_both_percentages_match.all()
    assert edges.left.str[0].isin(list("JKLP")).all() and edges.right.str[0].isin(list("JKLP")).all()
    assert edges.focal_delta.eq(0).all() and edges.all_duplicate_family_delta.eq(0).all()
    assert edges.identical_singleton_label_different_id.eq(0).all()
    affected = edges[edges.recovered_singleton_alias_pairs > 0]
    assert len(affected) == 9 and affected.recovered_singleton_alias_pairs.sum() == 65
    assert np.isclose((affected.recovered_singleton_alias_pairs / affected.paper_overlap).max(), 12/1245)
    ranked = table("fig1_literal_code_audit.tsv")
    exact = abs(ranked.reconstructed - ranked.published) < 1e-10
    assert len(ranked) == 1915 and exact.sum() == 1371
    assert exact[ranked.cluster.str.startswith("J_")].sum() == 116
    assert exact[ranked.cluster.str.startswith("K_")].sum() == 125
    fig7 = table("fig7_model_spot_axes.tsv")
    assert len(fig7) == 10 and fig7.raw_filtered_spots.sum() == 28965
    assert fig7.raw_state_order_identical.all() and fig7.raw_group_order_identical.all()
    assert fig7.group_columns.eq(5).all() and fig7.state_columns.eq(57).all()
    print("Verified original supplement, traced partial repair, Figure 1 and Figure 7 results.")


if __name__ == "__main__":
    main()
