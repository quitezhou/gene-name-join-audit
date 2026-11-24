"""v2 figures: Fig 1 (mechanism schematic + skin replay), Fig 2 (PBMC, prostate, adipose), Supp Fig S1."""
import json, sys
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Patch
from matplotlib.colors import BoundaryNorm, ListedColormap
import numpy as np, pandas as pd

V2 = Path(__file__).resolve().parent
HERE = V2.parents[1]
sys.path.insert(0, str(HERE))
import report_figures as rf

RED, GREEN, GREY, INK = "#b4462f", "#1f7a6d", "#8b959c", "#22303a"

def save(fig, name):
    fig.savefig(V2 / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(V2 / f"{name}.png", dpi=220, bbox_inches="tight")
    plt.close(fig)

def table(ax, x, y, title, rows, w=4.1, h=.5, highlight=None):
    ax.text(x, y + .12, title, fontsize=8.4, weight="bold", color=INK, va="bottom")
    out = []
    for k, (gid, name, keep) in enumerate(rows):
        yy = y - (k + 1) * h
        ax.add_patch(FancyBboxPatch((x, yy + .04), w, h * .82, boxstyle="round,pad=0.01,rounding_size=0.05",
                                    fc="#ffffff" if keep else "#f1f2f2", ec="#c3c9cd", lw=.6))
        col = INK if keep else "#a9b0b4"
        ax.text(x + .1, yy + h * .45, gid, fontsize=7.4, va="center", color=col, family="DejaVu Sans Mono")
        ax.text(x + w - .1, yy + h * .45, name if keep else "filtered out", fontsize=7.6, va="center", ha="right",
                color=col, style="normal" if keep else "italic",
                weight="bold" if (highlight and name == highlight and keep) else "normal")
        out.append(yy + h * .45)
    return out

def schematic(ax):
    ax.set_xlim(0, 14.2); ax.set_ylim(.1, 4.25); ax.axis("off")
    src = [("…271147", "ARMCX5-GPRASP2", True), ("…286237", "ARMCX5-GPRASP2", True), ("…", "other genes", True)]
    table(ax, 0.0, 3.85, "Shared source table", src)
    ax.text(0.0, 2.2, "Same Ensembl IDs, symbols\nand row order in both samples", fontsize=7.3, color=GREY, va="top")
    full = [("…271147", "ARMCX5-GPRASP2", True), ("…286237", "ARMCX5-GPRASP2-1", True)]
    sub = [("…271147", "", False), ("…286237", "ARMCX5-GPRASP2", True)]
    yf = table(ax, 5.0, 3.85, "A. Name both retained rows", full, highlight="ARMCX5-GPRASP2")
    ys = table(ax, 5.0, 1.75, "B. Name the retained row", sub, highlight="ARMCX5-GPRASP2")
    for y0 in (3.35, 1.25):
        ax.annotate("", xy=(4.85, y0), xytext=(4.2, 2.85), arrowprops=dict(arrowstyle="->", color=GREY, lw=.8))
    ax.annotate("", xy=(9.2, ys[1]), xytext=(9.2, yf[0]),
                arrowprops=dict(arrowstyle="<->", color=RED, lw=1.4, connectionstyle="arc3,rad=-0.45"))
    ax.text(10.0, 3.7, "Join by name", fontsize=8.4, weight="bold", color=RED, va="top")
    ax.text(10.0, 3.35, "ARMCX5-GPRASP2 pairs\n…271147 with …286237.\nEvery name matches.", fontsize=7.5, color=RED, va="top")
    ax.text(10.0, 1.75, "Join by Ensembl ID", fontsize=8.4, weight="bold", color=GREEN, va="top")
    ax.text(10.0, 1.4, "…286237 pairs with itself.", fontsize=7.5, color=GREEN, va="top")

def skin(bx):
    source = HERE.parent.parent
    order = json.loads((source / "tangram_skin_locus_adjudication.json").read_text())
    pairs = pd.DataFrame(order["pairs"])
    cases = pairs[pairs["reference"] != "POOLED_12_REFERENCES"]
    hits = cases[cases["mode"] == "filter_then_unique"].pivot(
        index="reference", columns="spatial", values="strict_training_conflicts").sort_index().sort_index(axis=1)
    assert hits.shape == (12, 16) and (hits > 0).to_numpy().sum() == 129
    assert (cases.loc[cases["mode"] != "filter_then_unique", "strict_training_conflicts"] == 0).all()
    colors = ["#eff2f1", "#91c7c4", "#3a8d8b", "#175252"]
    bx.imshow(hits, cmap=ListedColormap(colors), norm=BoundaryNorm([-.5,.5,1.5,2.5,3.5], 4),
              interpolation="nearest", aspect="auto")
    bx.set(xticks=[0, 3, 7, 11, 15], xticklabels=[1, 4, 8, 12, 16], xlabel="Spatial section",
           yticks=[0, 3, 7, 11], yticklabels=[1, 4, 8, 12], ylabel="Single-cell library")
    bx.set_xticks(np.arange(-.5, 16), minor=True); bx.set_yticks(np.arange(-.5, 12), minor=True)
    bx.grid(which="minor", color="white", lw=.65); bx.tick_params(which="minor", length=0)
    bx.set_title("b   Skin: filter, then name (counterfactual)", loc="left", fontsize=9.5)
    bx.legend(handles=[Patch(facecolor=c, edgecolor="#a9b9b9", linewidth=.4, label=str(i))
                       for i, c in enumerate(colors)], title="Disjoint-locus training keys per pair",
              ncol=4, loc="upper center", bbox_to_anchor=(.5, -.2), frameon=False, fontsize=7.5, title_fontsize=7.8)

def figure1():
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(10.2, 3.9), layout="constrained",
                                 gridspec_kw={"width_ratios": [1.7, 1]})
    schematic(ax)
    ax.set_title("a   Filtering before naming changes what a shared name denotes", loc="left", fontsize=9.5)
    skin(bx)
    save(fig, "figure_1_v2")

def prostate(ax):
    full = rf.load_result(rf.DATA / "all_results.tsv.gz")
    d = full[full.match == "one_release_label_same_id"].copy()
    d["status"] = "retained"
    d.loc[d.local_call & ~d.aligned_call, "status"] = "lost"
    d.loc[~d.local_call & d.aligned_call, "status"] = "gained"
    d.loc[d.local_call & d.aligned_call & (np.sign(d.local_log2FoldChange) !=
                                          np.sign(d.aligned_log2FoldChange)), "status"] = "reversed"
    n = d.status.value_counts().to_dict()
    for status, label, marker in [("retained", "Other", "o"), ("lost", "Call lost", "o"),
                                  ("gained", "Call gained", "^"), ("reversed", "Sign reversed", "D")]:
        sub = d[d.status == status].dropna(subset=["local_log2FoldChange", "aligned_log2FoldChange"])
        ax.scatter(sub.local_log2FoldChange, sub.aligned_log2FoldChange, s=13, c=rf.COLORS[status],
                   alpha=.8, label=f"{label} ({len(sub)})", marker=marker, linewidths=0)
    ax.axline((0, 0), slope=1, ls="--", lw=.8, c="#677079", zorder=0)
    ax.axhline(0, c="#bbc1c7", lw=.7, zorder=0)
    ax.set(xlabel="Joined by name: treatment log₂FC", ylabel="Joined by Ensembl ID: log₂FC")
    ax.set_title("b   Prostate: 251 renamed records", loc="left", fontsize=9.5)
    ax.legend(frameon=False, fontsize=7.2, loc="upper right", handletextpad=.2)
    return n

def figure2():
    fig = plt.figure(figsize=(10.2, 6.2), layout="constrained")
    gs = fig.add_gridspec(2, 2, height_ratios=[1.15, .8], width_ratios=[1.15, 1])
    ax, px, bx = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, :])
    # reuse PBMC + adipose drawing by rendering original function into our axes
    pbmc = pd.read_csv(HERE / "prospective_holdout/summary.tsv", sep="\t")
    within = pbmc[pbmc.comparison.isin(["ctrl", "stim"])]
    hits = within.pivot(index="threshold", columns="comparison", values="pairs_with_mismatch")
    locus = pd.read_csv(HERE / "prospective_holdout/locus_summary.tsv", sep="\t")
    sel = locus[(locus.release == "v19") & (locus.category == "disjoint") & locus.comparison.isin(hits.columns)
                & locus.threshold.astype(str).isin(hits.index.astype(str))]
    strict = sel.assign(threshold=sel.threshold.astype(int)).pivot(
        index="threshold", columns="comparison", values="donor_sample_pairs").reindex(hits.index).fillna(0).astype(int)
    x = np.arange(len(hits)); width = .34
    for field, off, color, label in [("ctrl", -.5, "#27696e", "Unstimulated"), ("stim", .5, "#bc7548", "IFN-β stimulated")]:
        counts = hits[field].to_numpy()
        bars = ax.bar(x + width*off, counts, width, color=color, label=label)
        ax.bar_label(bars, labels=[str(v) for v in counts], padding=2, fontsize=7.5)
        dist = strict[field].to_numpy()
        ax.bar(x + width*off, dist, width, hatch="////", edgecolor="white", linewidth=0, facecolor="none")
        for pos, v in zip(x + width*off, dist):
            if v:
                ax.text(pos, v/2, str(v), ha="center", va="center", color="white", fontsize=7.5, weight="bold",
                        bbox=dict(facecolor=color, edgecolor="none", pad=.5))
    ax.set(xticks=x, xticklabels=["No filter", "≥1", "≥3", "≥10"], ylim=(0, 31), yticks=[0, 7, 14, 21, 28],
           ylabel="Donor pairs with a\nsame-name, different-ID match (of 28)",
           xlabel="Minimum cells detecting the gene in each donor")
    ax.set_title("a   Blood: prespecified filtering thresholds", loc="left", fontsize=9.5)
    h, l = ax.get_legend_handles_labels()
    h.append(Patch(facecolor="white", edgecolor="#62686b", hatch="////", linewidth=.4))
    ax.legend(h, [*l, "Disjoint loci"], frameon=False, fontsize=7.5, loc="upper right")
    n = prostate(px)
    atlas = json.loads((HERE / "marker_resolution/contrast_bounds.json").read_text())["all_verified"]
    base, partial = atlas["released"]["contrast_pp"], atlas["partial"]["contrast_pp"]
    bx.axvline(0, color="#a8afb2", lw=.8, zorder=0)
    for y, label, pt, color in [(3, "Published row", base, "#79512f"), (2, "Source-verified cells realigned", partial, "#236b70")]:
        bx.scatter([pt], [y], c=color, s=28, zorder=2)
        bx.text(pt + .35, y, f"{pt:.2f}", va="center", ha="left", fontsize=7.5, color=color)
    for y, label, field in [(1, "Unresolved cells: either observed row", "two_row_cell_interval_pp"),
                            (0, "Unresolved cells: one row per source stratum", "two_row_source_stratum_interval_pp")]:
        lo, hi = atlas[field]
        bx.plot([lo, hi], [y, y], color="#424c54", lw=2.3, marker="|", ms=7)
        bx.text(hi + .35, y, f"{lo:.2f} to {hi:.2f}", va="center", ha="left", fontsize=7.5)
    bx.set(xlim=(-16, 20), ylim=(-.6, 3.6), yticks=[3, 2, 1, 0], xticks=range(-15, 20, 5),
           yticklabels=["Published row", "Source-verified cells realigned",
                        "Unresolved cells: either observed row", "Unresolved cells: one row per source stratum"],
           xlabel="POLR2J3.1 detection, T-cell cluster 3 minus other T cells (percentage points)")
    bx.spines["left"].set_visible(False); bx.tick_params(axis="y", length=0)
    bx.set_title("c   Adipose atlas: a published marker traced to source records", loc="left", fontsize=9.5)
    save(fig, "figure_2_v2")
    return n

def supp_s1():
    import types
    # Draw only the cross-table panel from the original function by re-running it and keeping panel a.
    source = HERE.parent.parent
    panel = json.loads((source / "feature_panel_result.json").read_text())
    axes = panel["dataset_summaries"]
    ex = ["PBMC3K_hg19", "PBMC8K_GRCh38_1.2.0", "PBMC10K_GRCh38_3.0.0",
          "Prostate_spatial_NEADT12_GRCh38_2020_A", "Adrenal_normal_spatial_NAd4"]
    disp = ["A   hg19", "B   GRCh38 1.2", "C   GRCh38 3.0", "D   GRCh38 2020-A", "E   FFPE panel"]
    groups = {axes[k]["ordered_id_symbol_sha256"]: i for i, k in enumerate(ex)}
    counts = np.full((5, 5), np.nan)
    for p in panel["pairwise_results"]:
        i, j = sorted(groups[axes[k]["ordered_id_symbol_sha256"]] for k in (p["first"], p["second"]))
        if i != j: counts[i, j] = p["duplicate_driven_same_name_different_id"]
    fig, ax = plt.subplots(figsize=(4.6, 3.8), layout="constrained")
    ax.imshow(counts, cmap="Blues", vmin=0, vmax=43)
    for i, j in zip(*np.where(np.isfinite(counts))):
        v = int(counts[i, j]); ax.text(j, i, str(v), ha="center", va="center", color="white" if v >= 27 else "#143547")
    ax.set(yticks=range(5), yticklabels=disp, xticks=range(5), xticklabels=list("ABCDE"), xlabel="Ordered feature table")
    save(fig, "supp_figure_S1")

if __name__ == "__main__":
    plt.rcParams.update({"font.size": 8.5, "font.family": "DejaVu Sans", "pdf.fonttype": 42,
                         "axes.spines.top": False, "axes.spines.right": False})
    figure1(); print(figure2()); supp_s1()
