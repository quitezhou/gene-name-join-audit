"""Generate the two Brief Report figures and reconcile reported denominators."""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
CASE = HERE.parent / "sentinel_case"
HWAT = HERE.parent / "hwat_verification"
DATA = HERE / "full_S11"
COLORS = dict(retained="#8896a3", lost="#b95639", reversed="#7f3089", gained="#197f75")


def load_result(filename):
    return pd.read_csv(filename, sep="\t", index_col="gene")


def summarize(full):
    groups = {}
    for kind, frame in list(full.groupby("match")) + [("all_rows", full)]:
        local, aligned, paper = (frame[k + "_call"] for k in ("local", "aligned", "published"))
        reversed_sign = np.sign(frame.local_log2FoldChange) != np.sign(frame.aligned_log2FoldChange)
        rho = frame.local_log2FoldChange.corr(frame.published_log2FoldChange, method="spearman")
        groups[kind] = {
            "rows": len(frame), "changed_count_vectors": int((frame.changed_count_samples > 0).sum()),
            "published_calls": int(paper.sum()), "local_calls": int(local.sum()),
            "aligned_calls": int(aligned.sum()), "lost": int((local & ~aligned).sum()),
            "gained": int((~local & aligned).sum()),
            "published_and_local_but_lost": int((paper & local & ~aligned).sum()),
            "opposite_sign_passing_both": int((local & aligned & reversed_sign).sum()),
            "spearman_published_vs_local_lfc": rho if pd.notna(rho) else None}
    return groups


def save(fig, name):
    fig.savefig(HERE / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(HERE / f"{name}.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_mechanism():
    from matplotlib.colors import BoundaryNorm, ListedColormap
    from matplotlib.patches import Patch

    source = HERE.parent.parent
    panel = json.loads((source / "feature_panel_result.json").read_text())
    axes = panel["dataset_summaries"]
    examples = ["PBMC3K_hg19", "PBMC8K_GRCh38_1.2.0", "PBMC10K_GRCh38_3.0.0",
                "Prostate_spatial_NEADT12_GRCh38_2020_A", "Adrenal_normal_spatial_NAd4"]
    displays = ["A   hg19", "B   GRCh38 1.2", "C   GRCh38 3.0", "D   GRCh38 2020-A", "E   FFPE panel"]
    groups = {axes[k]["ordered_id_symbol_sha256"]: i for i, k in enumerate(examples)}
    assert len(groups) == len({r["ordered_id_symbol_sha256"] for r in axes.values()}) == 5
    counts = np.full((5, 5), np.nan)
    for pair in panel["pairwise_results"]:
        i, j = sorted(groups[axes[k]["ordered_id_symbol_sha256"]] for k in (pair["first"], pair["second"]))
        n = pair["duplicate_driven_same_name_different_id"]
        if i == j:
            assert n == 0
        else:
            assert np.isnan(counts[i, j]) or counts[i, j] == n
            counts[i, j] = n
    assert np.isfinite(counts).sum() == 10
    order = json.loads((source / "tangram_skin_locus_adjudication.json").read_text())
    pairs = pd.DataFrame(order["pairs"])
    cases = pairs[pairs["reference"] != "POOLED_12_REFERENCES"]
    hits = cases[cases["mode"] == "filter_then_unique"].pivot(
        index="reference", columns="spatial", values="strict_training_conflicts").sort_index().sort_index(axis=1)
    assert hits.shape == (12, 16) and hits.notna().all().all() and (hits > 0).to_numpy().sum() == 129
    assert (cases.loc[cases["mode"] != "filter_then_unique", "strict_training_conflicts"] == 0).all()

    fig, (ax, bx) = plt.subplots(1, 2, figsize=(9.4, 3.8), layout="constrained",
                                 gridspec_kw={"width_ratios": [1, 1.25]})
    ax.imshow(counts, cmap="Blues", vmin=0, vmax=43, interpolation="nearest")
    for i, j in zip(*np.where(np.isfinite(counts))):
        n = int(counts[i, j])
        ax.text(j, i, str(n), ha="center", va="center", color="white" if n >= 27 else "#143547", fontsize=10)
    ax.set(yticks=range(5), yticklabels=displays, xticks=range(5), xticklabels=list("ABCDE"),
           xlabel="Ordered feature table")
    ax.set_title("a   Reassigned names across feature tables", loc="left", fontsize=10.5)
    colors = ["#eff2f1", "#91c7c4", "#3a8d8b", "#175252"]
    bx.imshow(hits, cmap=ListedColormap(colors), norm=BoundaryNorm([-.5,.5,1.5,2.5,3.5], 4),
              interpolation="nearest", aspect="auto")
    bx.set(xticks=[0, 3, 7, 11, 15], xticklabels=[1, 4, 8, 12, 16], xlabel="Spatial section",
           yticks=[0, 3, 7, 11], yticklabels=[1, 4, 8, 12], ylabel="Single-cell reference")
    bx.set_xticks(np.arange(-.5, 16), minor=True)
    bx.set_yticks(np.arange(-.5, 12), minor=True)
    bx.grid(which="minor", color="white", lw=.65)
    bx.tick_params(which="minor", length=0)
    bx.set_title("b   Different loci after filtering before naming", loc="left", fontsize=10.5)
    bx.legend(handles=[Patch(facecolor=c, edgecolor="#a9b9b9", linewidth=.4, label=str(i))
                       for i, c in enumerate(colors)], title="Mismatched training keys per pair",
              ncol=4, loc="upper center", bbox_to_anchor=(.5, -.18), frameon=False, fontsize=8, title_fontsize=8.5)
    save(fig, "figure_1_mechanism")


def plot_prostate(reconciled):
    d = reconciled[reconciled.match == "one_release_label_same_id"].copy()
    d["status"] = "retained"
    d.loc[d.local_call & ~d.aligned_call, "status"] = "lost"
    d.loc[~d.local_call & d.aligned_call, "status"] = "gained"
    d.loc[d.local_call & d.aligned_call & (np.sign(d.local_log2FoldChange) !=
                                          np.sign(d.aligned_log2FoldChange)), "status"] = "reversed"
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(9, 4.1), layout="constrained",
                                 gridspec_kw={"width_ratios": [1.17, 1]})
    for status, label, marker in [("retained", "Other", "o"), ("lost", "Call lost", "o"),
                                  ("gained", "Call gained", "^"), ("reversed", "Direction reversed", "D")]:
        sub = d[d.status == status].dropna(subset=["local_log2FoldChange", "aligned_log2FoldChange"])
        ax.scatter(sub.local_log2FoldChange, sub.aligned_log2FoldChange, s=19, c=COLORS[status],
                   alpha=.8, label=f"{label} ({len(sub)})", marker=marker, linewidths=0)
    ax.axline((0, 0), slope=1, ls="--", lw=.8, c="#677079", zorder=0)
    ax.axhline(0, c="#bbc1c7", lw=.7, zorder=0)
    ax.set(xlabel="Local name · treatment log₂ fold change", ylabel="Shared ID · treatment log₂ fold change")
    ax.set_title("a   Renamed records in the test panel", loc="left", fontsize=11)
    ax.legend(frameon=False, fontsize=8, loc="upper center", bbox_to_anchor=(.50, -.19), ncol=2)
    # First alphabetical renamed record. Each point is one observed section.
    assert d.sort_index().index[0] == "ABITRAM"
    archive = np.load(CASE / "pseudobulk_counts.npz")
    strata = {"Untreated\nolder": [], "Treated\nolder": [], "Treated\nnewer": []}
    for rec in json.loads((CASE / "matrix_manifest.json").read_text()):
        acc, sample = rec["accession"], rec["sample"]
        source = CASE.parent.parent / "panel_sources/all_features" / f"{acc}_{sample}_features.tsv.gz"
        features = pd.read_csv(source, sep="\t", names=["id", "symbol", "type"])
        anchor = d.at["ABITRAM", "anchor"]
        ix = np.flatnonzero(features.id.eq(anchor))
        assert len(ix) == 1
        i = ix[0]
        measured = archive[acc + "_qc"][i] * 1e6 / rec["qc_UMI"]
        local = measured if features.symbol.iloc[i] == "ABITRAM" else 0.0
        if len(features) == 36601:
            key = "Treated\nnewer"
        else:
            key = "Untreated\nolder" if sample.startswith("TRNA_") else "Treated\nolder"
        strata[key].append((local, measured))
    colors = ["#9a543b", "#206f80"]
    for j, (label, rows) in enumerate(strata.items()):
        vals = np.asarray(rows)
        for i, col in enumerate(colors):
            x = j + (i - .5) * .27
            jitter = np.linspace(-.055, .055, len(rows))
            bx.scatter(x + jitter, vals[:, i], c=col, alpha=.6, s=12,
                       label=["Local ABITRAM", "Shared Ensembl ID"][i] if j == 0 else None)
            bx.plot([x-.08, x+.08], [np.median(vals[:, i])]*2, color=col, lw=2)
    bx.set(xticks=np.arange(3), xticklabels=list(strata), ylabel="ABITRAM counts per million")
    bx.set_title("b   Counts restored under FAM206A", loc="left", fontsize=11)
    bx.legend(frameon=False, fontsize=8, loc="upper center", bbox_to_anchor=(.50, -.19), ncol=2)
    save(fig, "figure_1_prostate")


def plot_atlas():
    df = pd.read_csv(HWAT / "adipocyte_detection_contrasts.tsv", sep="\t").set_index("symbol")
    marker = pd.read_csv(HERE / "hwat_model/tcell_published_marker_detection.tsv", sep="\t").set_index("arm")
    gene = df.loc["ARMCX5-GPRASP2"]
    published = marker.loc["released_object_printed_key"]
    repaired = marker.loc["released_object_source_verified_partial_id_repair"]
    assert (published[["in_n", "out_n"]].to_numpy() == [1154, 9240]).all()
    fig, axs = plt.subplots(1, 2, figsize=(8.4, 3.85), sharey=True, layout="constrained")
    panels = [
        ("a   ARMCX5-GPRASP2 in adipocytes", ["Library 1\n3,571 cells", "Library 4\n306 cells"],
         [("Bare name", gene[["bare_pct_1", "bare_pct_2"]]),
          ("Shared ID", gene[["shared_id_pct_1", "shared_id_pct_2"]])]),
        ("b   POLR2J3.1 in atlas T cells", ["Cluster 3\n1,154 cells", "Other T cells\n9,240 cells"],
         [("Published row", 100 * published[["in_fraction", "out_fraction"]]),
          ("Partial ID alignment", 100 * repaired[["in_fraction", "out_fraction"]])])]
    x = np.arange(2)
    for (title, labels, arms), ax in zip(panels, axs):
        for (key, heights), color, dx in zip(arms, ["#9a543b", "#206f80"], [-.16, .16]):
            ax.bar(x+dx, heights, width=.29, color=color, label=key)
            for xx, h in zip(x+dx, heights):
                ax.text(xx, h+.4, f"{h:.2f}", ha="center", fontsize=8)
        ax.set(xticks=x, xticklabels=labels, ylim=(0, 18))
        ax.set_title(title, loc="left", fontsize=11)
        ax.legend(frameon=False, fontsize=8, loc="upper center", bbox_to_anchor=(.5, -.19), ncol=2)
    axs[0].set_ylabel("Cells with detected signal (%)")
    save(fig, "figure_2_atlas")


def plot_replication():
    from matplotlib.patches import Patch

    pbmc = pd.read_csv(HERE / "prospective_holdout/summary.tsv", sep="\t")
    within = pbmc[pbmc.comparison.isin(["ctrl", "stim"])]
    hits = within.pivot(index="threshold", columns="comparison", values="pairs_with_mismatch")
    assert hits.index.tolist() == [0, 1, 3, 10] and within.pairs.eq(28).all()
    locus = pd.read_csv(HERE / "prospective_holdout/locus_summary.tsv", sep="\t")
    selected = locus[(locus.release == "v19") & (locus.category == "disjoint")
                     & locus.comparison.isin(hits.columns) & locus.threshold.astype(str).isin(hits.index.astype(str))]
    strict = selected.assign(threshold=selected.threshold.astype(int)).pivot(
        index="threshold", columns="comparison", values="donor_sample_pairs").loc[hits.index]
    atlas = json.loads((HERE / "marker_resolution/contrast_bounds.json").read_text())["all_verified"]
    base = atlas["released"]["contrast_pp"]
    partial = atlas["partial"]["contrast_pp"]
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(9.2, 4.05), layout="constrained",
                                 gridspec_kw={"width_ratios":[1.05, 1]})
    x = np.arange(len(hits))
    width = .34
    for field, offset, color, label in [("ctrl", -.5, "#27696e", "Unstimulated"),
                                        ("stim", .5, "#bc7548", "IFN-β stimulated")]:
        counts = hits[field].to_numpy()
        bars = ax.bar(x+width*offset, counts, width, color=color, label=label)
        ax.bar_label(bars, labels=[str(n) for n in counts], padding=3, fontsize=8)
        distinct = strict[field].to_numpy()
        assert (distinct <= counts).all()
        ax.bar(x+width*offset, distinct, width, hatch="////", edgecolor="white",
               linewidth=0, facecolor="none")
        for pos, value in zip(x+width*offset, distinct):
            if value:
                ax.text(pos, value/2, str(value), ha="center", va="center", color="white",
                        fontsize=8, weight="bold", bbox=dict(facecolor=color, edgecolor="none", pad=.5))
    ax.set(xticks=x, xticklabels=["No filter", "≥1", "≥3", "≥10"], ylim=(0, 31),
           yticks=[0, 7, 14, 21, 28], ylabel="Donor pairs with a different-ID match (of 28)",
           xlabel="Minimum cells detecting the gene in each donor")
    ax.set_title("a   PBMC gene filtering in eight donors", loc="left", fontsize=10.5)
    handles, labels = ax.get_legend_handles_labels()
    handles.append(Patch(facecolor="white", edgecolor="#62686b", hatch="////", linewidth=.4))
    ax.legend(handles, [*labels, "Disjoint loci"], frameon=False, fontsize=8, loc="upper right")
    bx.axvline(0, color="#a8afb2", lw=.8, zorder=0)
    # Labels are printed above each mark to keep the panel readable in one column.
    for y, label, point, color in [(3, "Published row", base, "#79512f"),
                                   (2, "Partial ID repair", partial, "#236b70")]:
        bx.scatter([point], [y], c=color, s=32, zorder=2)
        bx.text(-15, y+.09, label, va="bottom", fontsize=9)
        bx.text(point, y-.07, f"{point:.2f}", va="top", ha="center", fontsize=8, color=color)
    for y,label,field in [(1,"Either observed row per cell","two_row_cell_interval_pp"),
                         (0,"One row per source stratum","two_row_source_stratum_interval_pp")]:
        lo, hi = atlas[field]
        bx.plot([lo,hi],[y,y], color="#424c54", lw=2.5, marker="|", ms=7)
        bx.text(-15, y+.09, label, va="bottom", fontsize=9)
        for value in (lo,hi):
            bx.text(value, y-.08, f"{value:.2f}", ha="center", va="top", fontsize=8)
    bx.set(xlim=(-16,16), ylim=(-.43,3.5), yticks=[], xlabel="Cluster 3 minus other T cells (percentage points)")
    bx.spines["left"].set_visible(False)
    bx.set_title("b   Adipose T-cell marker detection", loc="left", fontsize=10.5)
    save(fig, "figure_2_replication")


if __name__ == "__main__":
    plt.rcParams.update({"font.size": 9, "font.family": "DejaVu Sans", "pdf.fonttype": 42,
                         "axes.spines.top": False, "axes.spines.right": False})
    full = load_result(DATA / "all_results.tsv.gz")
    common = load_result(DATA / "common_axis/all_results.tsv.gz")
    old = load_result(DATA / "common_axis/old_release/all_results.tsv.gz")
    summary = {key: summarize(x) for key, x in [("main", full), ("common", common), ("old", old)]}
    (DATA / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False))
    plot_mechanism()
    plot_prostate(full)
    plot_atlas()
    plot_replication()
