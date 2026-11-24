"""Exhaustively reconcile S11 labels against two deposited feature axes.

Run in the pinned sentinel_case PyDESeq2 environment. --inventory skips fitting.
"""
import argparse
import json
from importlib.metadata import version
from pathlib import Path

import numpy as np
import pandas as pd
from anndata.utils import make_index_unique

HERE = Path(__file__).resolve().parent
CASE = HERE.parents[1] / "sentinel_case"
PROJECT = CASE.parent.parent


def inputs():
    paper = pd.read_excel(CASE / "source/MOESM13.xlsx", index_col="gene")
    manifest = json.loads((CASE / "matrix_manifest.json").read_text())
    existing = json.loads((PROJECT / "prostate_spatial_versions_result.json").read_text())
    switches = {row["local_unique_name"]: row for row in
                existing["comparison_after_per_table_anndata_uniquification"]["mismatch_examples"]}
    archive = np.load(CASE / "pseudobulk_counts.npz")
    axes, samples = {}, []
    for rec in manifest:
        gene_file = PROJECT / "panel_sources/all_features" / f"{rec['accession']}_{rec['sample']}_features.tsv.gz"
        features = pd.read_csv(gene_file, sep="\t", names=["id", "symbol", "type"])
        assert features.id.is_unique and features.type.eq("Gene Expression").all()
        local = make_index_unique(pd.Index(features.symbol))
        axis = pd.DataFrame({"id": features.id.values, "symbol": features.symbol.values}, index=local)
        release = {33538: "old", 36601: "new"}[len(axis)]
        if release in axes:
            assert axes[release].equals(axis)
        else:
            axes[release] = axis
        counts = archive[rec["accession"] + "_qc"]
        assert counts.size == len(axis) and counts.sum() == rec["qc_UMI"]
        samples.append((rec["sample"], release, counts))
    assert len(samples) == 39 and len(paper) == 17696 and paper.index.is_unique
    return paper, axes, switches, samples


def adjudicate(labels, axes, switches):
    ids = {k: set(axis.id) for k, axis in axes.items()}
    names = {k: dict(zip(axis.id, axis.index)) for k, axis in axes.items()}
    rows = []
    for label in labels:
        old, new = (axes[k].id.get(label) for k in ("old", "new"))
        anchor = None
        if old is None and new is None:
            kind = "label_missing_both"
        elif old == new:
            kind, anchor = "identical", old
        elif old is None or new is None:
            candidate = old or new
            if candidate in ids["old"] and candidate in ids["new"]:
                kind, anchor = "one_release_label_same_id", candidate
            else:
                kind = "id_missing_in_one_release"
        elif label in switches:
            row = switches[label]
            common = set(row["first_symbol_group"]) & set(row["second_symbol_group"])
            assert len(common) == 1 and common <= ids["old"] & ids["new"]
            kind, anchor = "shared_label_duplicate_switch", common.pop()
        else:
            kind = "different_ids_no_unique_anchor"
        rows.append(dict(gene=label, old_local_id=old, new_local_id=new, match=kind,
                         anchor=anchor, anchor_old_label=names["old"].get(anchor),
                         anchor_new_label=names["new"].get(anchor)))
    lookup = pd.DataFrame(rows).set_index("gene")
    lookup["nonunique_panel_anchor"] = lookup.anchor.notna() & lookup.anchor.duplicated(keep=False)
    lookup["substitute"] = lookup.match.isin(["one_release_label_same_id", "shared_label_duplicate_switch"]) & ~lookup.nonunique_panel_anchor
    return lookup


def make_counts(lookup, axes, samples):
    baseline, repaired = {}, {}
    replace = lookup.index[lookup.substitute]
    for name, release, vector in samples:
        axis = axes[release]
        values = pd.Series(vector, index=axis.index)
        baseline[name] = values.reindex(lookup.index, fill_value=0)
        repaired[name] = baseline[name].copy()
        repaired[name].loc[replace] = pd.Series(vector, index=axis.id).loc[lookup.loc[replace, "anchor"]].to_numpy()
    return pd.DataFrame(baseline).T, pd.DataFrame(repaired).T


def fit(counts, factors, n_cpus):
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    from pydeseq2.default_inference import DefaultInference

    metadata = pd.DataFrame({"treatment": [s.split("_")[0] for s in counts.index]}, index=counts.index)
    inference = DefaultInference(n_cpus=n_cpus)
    dds = DeseqDataSet(counts=counts, metadata=metadata, design_factors="treatment",
                       ref_level=["treatment", "TRNA"], inference=inference)
    dds.obsm["size_factors"] = factors.copy()
    dds.layers["normed_counts"] = dds.X / factors[:, None]
    dds.fit_genewise_dispersions()
    dds.fit_dispersion_trend()
    dds.fit_dispersion_prior()
    dds.fit_MAP_dispersions()
    dds.fit_LFC()
    dds.calculate_cooks()
    if dds.refit_cooks:
        dds.refit()
    stats = DeseqStats(dds, inference=inference, contrast=["treatment", "NEADT", "TRNA"])
    stats.summary()
    return stats.results_df


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inventory", action="store_true")
    parser.add_argument("--common-axis", action="store_true",
                        help="Post hoc sensitivity: test only single, shared-ID anchors")
    parser.add_argument("--old-release", action="store_true",
                        help="Post hoc sensitivity: keep the 28 samples with cross-treatment overlap")
    parser.add_argument("--cores", type=int, default=4)
    args = parser.parse_args()
    paper, axes, switches, samples = inputs()
    lookup = adjudicate(paper.index, axes, switches)
    baseline, repaired = make_counts(lookup, axes, samples)
    outdir = HERE
    if args.common_axis:
        keep = lookup.anchor.notna() & ~lookup.nonunique_panel_anchor
        lookup = lookup.loc[keep].copy()
        baseline = baseline.loc[:, keep]
        repaired = repaired.loc[:, keep]
        outdir = outdir / "common_axis"
    if args.old_release:
        names = [sample for sample, release, _ in samples if release == "old"]
        baseline, repaired = baseline.loc[names], repaired.loc[names]
        outdir = outdir / "old_release"
    outdir.mkdir(exist_ok=True, parents=True)
    lookup["changed_count_samples"] = (baseline != repaired).sum()
    assert baseline.loc[:, ~lookup.substitute].equals(repaired.loc[:, ~lookup.substitute])
    for arm, counts in (("local", baseline), ("aligned", repaired)):
        for group in ("TRNA", "NEADT"):
            lookup[f"{arm}_{group}_UMI"] = counts.loc[counts.index.str.startswith(group + "_")].sum()
    lookup.to_csv(outdir / "all_feature_keys.tsv.gz", sep="\t")
    print(lookup.groupby(["match", "nonunique_panel_anchor", "substitute"]).size(), flush=True)
    if args.inventory:
        return
    assert version("pydeseq2") == "0.4.4"
    from pydeseq2.preprocessing import deseq2_norm_fit, deseq2_norm_transform
    background_cols = lookup.index[~lookup.substitute]
    background = baseline[background_cols].to_numpy()
    logmeans, keep = deseq2_norm_fit(background)
    _, factors = deseq2_norm_transform(background, logmeans, keep)
    pd.Series(factors, index=baseline.index, name="common_size_factor").to_csv(outdir / "size_factors.tsv", sep="\t")
    out = lookup.copy()
    for arm, table in (("published", paper),
                       ("local", fit(baseline, factors, args.cores)),
                       ("aligned", fit(repaired, factors, args.cores))):
        for key in ("log2FoldChange", "lfcSE", "padj"):
            out[f"{arm}_{key}"] = table[key]
        out[f"{arm}_call"] = (table.log2FoldChange.abs() >= 1) & (table.padj < 0.05)
    out.to_csv(outdir / "all_results.tsv.gz", sep="\t")
    out.loc[out.substitute].to_csv(outdir / "target_results.tsv", sep="\t")


if __name__ == "__main__":
    main()
