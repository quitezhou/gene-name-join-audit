"""Verify the enclosed tables and render the two figures without network access."""
import csv
import hashlib
import gzip
import json
import os
import re
import runpy
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency

ROOT = Path(__file__).resolve().parent
A = ROOT / "audit_20261123"
F = A / "gb_followup"


def local_only(event, args):
    if event in {"socket.connect", "socket.getaddrinfo", "socket.sendto"}:
        raise RuntimeError("The offline review attempted a network operation")


def verify_payload():
    for row in csv.DictReader((ROOT / "MANIFEST.tsv").open(), delimiter="\t"):
        contents = (ROOT / row["path"]).read_bytes()
        assert len(contents) == int(row["bytes"]), row["path"]
        assert hashlib.sha256(contents).hexdigest() == row["sha256"], row["path"]
    dest = ROOT / "panel_sources/all_features"
    dest.mkdir(exist_ok=True)
    for row in json.loads((ROOT / "axes/index.json").read_text()):
        data = (ROOT / "axes" / (row["blob"] + ".tsv.gz")).read_bytes()
        reconstructed = bytes.fromhex(row["header"]) + data[10:]
        assert hashlib.sha256(reconstructed).hexdigest() == row["sha256"]
        (dest / row["name"]).write_bytes(reconstructed)
    methods = (F / "Supplementary_Methods.md").read_text(encoding="utf8")
    for link in re.findall(r"\]\(([^)]+)\)", methods):
        if "://" not in link:
            assert (F / link).is_file(), link
    initial = (F / "cohort_followup.tsv").read_bytes()
    runpy.run_path(str(F / "update_cohort.py"), run_name="__main__")
    assert (F / "cohort_followup.tsv").read_bytes() == initial
    for row in csv.DictReader((F / "cohort_followup.tsv").open(encoding="utf8"), delimiter="\t"):
        if row["evidence_folder"]:
            assert (F / row["evidence_folder"]).is_dir()


def table(path):
    return pd.read_csv(path, sep="\t")


def verify_mechanism():
    # Run both entry points inside this process to keep our network guard active.
    folder = F / "systematic_mechanism"
    old_path = list(sys.path)
    old_args = list(sys.argv)
    try:
        sys.path.insert(0, str(folder))
        with tempfile.TemporaryDirectory() as scratch:
            script = folder / "check_panel.py"
            sys.argv = [str(script), "--output-dir", scratch]
            runpy.run_path(str(script), run_name="__main__")
            for name in ["panel_pairs.tsv", "panel_events.tsv.gz", "panel_summary.json"]:
                expected, actual = folder / name, Path(scratch) / name
                if name.endswith(".gz"):
                    assert gzip.decompress(actual.read_bytes()) == gzip.decompress(expected.read_bytes())
                else:
                    assert actual.read_text(encoding="utf-8") == expected.read_text(encoding="utf-8")
        script = folder / "verify_order.py"
        sys.argv = [str(script)]
        checks = runpy.run_path(str(script), run_name="__main__")
        verify_mechanism_plot_inputs(folder, checks["distinct_loci"])
        with tempfile.TemporaryDirectory() as scratch:
            # The lung replay invokes the public checker in two new interpreters.
            # Apply the same network guard in these child processes at startup.
            hook = Path(scratch) / "sitecustomize.py"
            hook.write_text(
                "import sys\n"
                "def local_only(event, args):\n"
                "    if event in {'socket.connect', 'socket.getaddrinfo', 'socket.sendto'}:\n"
                "        raise RuntimeError('Network access during offline review')\n"
                "sys.addaudithook(local_only)\n", encoding="utf-8")
            previous = os.environ.get("PYTHONPATH")
            os.environ["PYTHONPATH"] = scratch + (os.pathsep + previous if previous else "")
            try:
                script = folder / "check_lung_join.py"
                out = Path(scratch) / "lung_join_check.json"
                sys.argv = [str(script), "--output", str(out)]
                runpy.run_path(str(script), run_name="__main__")
                assert json.loads(out.read_text()) == json.loads((folder / out.name).read_text())
            finally:
                if previous is None:
                    del os.environ["PYTHONPATH"]
                else:
                    os.environ["PYTHONPATH"] = previous
    finally:
        sys.path[:] = old_path
        sys.argv[:] = old_args


def verify_mechanism_plot_inputs(folder, distinct_loci):
    # Verify Figure 1 cell by cell against the independent axis and locus audits.
    panel = json.loads((ROOT / "feature_panel_result.json").read_text())
    axes = panel["dataset_summaries"]
    observed = {frozenset((r.left_axis, r.right_axis)): r.different_ids_duplicate
                for r in table(folder / "panel_pairs.tsv").itertuples()}
    for pair in panel["pairwise_results"]:
        key = frozenset(axes[pair[side]]["ordered_id_symbol_sha256"] for side in ("first", "second"))
        assert pair["duplicate_driven_same_name_different_id"] == (observed[key] if len(key) == 2 else 0)
    with gzip.open(ROOT / "gencode.v32.gene_metadata.json.gz", "rt", encoding="utf-8") as fh:
        genes = json.load(fh)["gene_records"]
    counts = json.loads((ROOT / "tangram_skin_real_counts_result.json").read_text())
    plotted = json.loads((ROOT / "tangram_skin_locus_adjudication.json").read_text())
    def index(row):
        return row["reference"], row["spatial"], row["mode"]
    expected = {index(row): sorted(row["strict_training_keys"]) for row in plotted["pairs"]}
    assert len(expected) == len(counts["pairs"]) == 624
    for row in counts["pairs"]:
        hits = [c["key"] for c in row["conflicts"] if c["enters_default_training"]
                and distinct_loci(genes, c["reference_id"], c["spatial_id"])]
        assert sorted(hits) == expected.pop(index(row))
    assert not expected


def prostate():
    p = F / "full_S11"
    all_results = table(p / "all_results.tsv.gz").set_index("gene")
    m = all_results[all_results.match.eq("one_release_label_same_id")]
    assert len(all_results) == 17696 and len(m) == 251 and (m.changed_count_samples > 0).sum() == 245
    for arm in ["local", "aligned", "published"]:
        above = m[f"{arm}_log2FoldChange"].abs().ge(1) & m[f"{arm}_padj"].lt(.05)
        assert np.array_equal(above, m[f"{arm}_call"])
    loss = m.local_call & ~m.aligned_call
    recovered = m.index[loss & m.published_call]
    opposite = m.index[m.local_call & m.aligned_call & (m.local_log2FoldChange * m.aligned_log2FoldChange < 0)]
    assert (m.local_call.sum(), m.aligned_call.sum(), loss.sum()) == (192, 29, 169)
    assert len(recovered) == 167 and (~m.local_call & m.aligned_call).sum() == 6
    assert m.loc[recovered, "local_TRNA_UMI"].eq(0).all() and m.loc[recovered, "aligned_TRNA_UMI"].gt(0).all()
    assert len(opposite) == 7 and m.loc[opposite, "published_call"].all()
    common = table(p / "common_axis/all_results.tsv.gz").set_index("gene")
    old = table(p / "common_axis/old_release/all_results.tsv.gz").set_index("gene")
    assert len(common) == 17495 and m.aligned_call.equals(common.loc[m.index, "aligned_call"])
    assert old.loc[opposite, "aligned_log2FoldChange"].lt(0).all() and old.loc[opposite, "aligned_call"].sum() == 5
    assert m.at["ABITRAM","aligned_TRNA_UMI"] == 5044
    marks = table(A / "sentinel_case/fig5i_marker_check.tsv")
    # Verify Figure 5i stars in the complete paired fits, not only the earlier targeted table.
    def stars(padj):
        return "***" if padj < .001 else "**" if padj < .01 else "*" if padj < .05 else "none"
    both = all_results.loc[marks.gene, ["local_padj", "aligned_padj"]]
    assert len(both) == 18 and (both.local_padj.map(stars) == both.aligned_padj.map(stars)).all()
    manifests = json.loads((A / "sentinel_case/matrix_manifest.json").read_text())
    assert len(manifests) == 39 and sum(m["qc_spots"] for m in manifests) == 105833


def atlas():
    p = F / "hwat_model"
    s = table(p / "tcell_published_marker_detection.tsv").set_index("arm")
    orig, fix = s.loc["released_object_printed_key"], s.loc["released_object_source_verified_partial_id_repair"]
    assert (orig.in_positive, orig.in_n, orig.out_positive, orig.out_n) == (151, 1154, 233, 9240)
    assert (fix.in_positive, fix.out_positive) == (163, 487)
    trace = table(p / "tcell_published_marker_id_trace.tsv")
    trace = trace[trace.include]
    assert len(trace) == 4 and trace.barcode_pairs.sum() == 1375
    for arm in ["in", "out"]:
        assert (trace["id_positive_" + arm] - trace["printed_key_positive_" + arm]).sum() == fix[arm + "_positive"] - orig[arm + "_positive"]
    matched = table(p / "fig1b_tested_edges.tsv")
    assert len(matched) == 277 and matched.weight_and_both_percentages_match.all()
    assert matched.all_duplicate_family_delta.eq(0).all()
    assert matched.recovered_singleton_alias_pairs.gt(0).sum() == 9
    adipo = table(A / "hwat_verification/adipocyte_detection_contrasts.tsv").set_index("symbol").loc["ARMCX5-GPRASP2"]
    assert (adipo.bare_positive_1, adipo.bare_positive_2, adipo.shared_id_positive_2) == (41, 47, 3)
    endpoints = table(F / "atlas_endpoint/paired_summary.tsv").set_index("model")
    assert endpoints.cells.eq(7899).all() and endpoints.loc[["Group","High","Low"], "top_changed"].tolist() == [3,3,8]
    official = table(F / "atlas_endpoint/official_api_verification.tsv")
    assert official.label_disagreements.eq(0).all()


def lung():
    folder = F / "fifth_case"
    raw = table(folder / "input_axes.tsv")
    assert raw.accession_tail.value_counts().to_dict() == {526:81, 530:40}
    summary = json.loads((folder / "summary.json").read_text())
    assert (summary["matched_local"], summary["different_id"], summary["same_id"]) == (24769, 221, 24548)
    discord = table(folder / "discordant_model_keys.tsv")
    dupe = discord[discord.cause == "duplicate_associated"]
    assert len(dupe) == 6 and len(discord) == 221
    assert ((discord.cause == "singleton_name") & (discord.gencode_v32_interval == "disjoint")).sum() == 140
    umi = ["Background_UMI", "Tumour_UMI", "Healthy_UMI"]
    assert dupe[umi].sum().sum() == 45797 and dupe.set_index("input_key").loc["ATXN7", umi].sum() == 3924
    idx = ["Patient_number", "Tissue", "pair"]
    counts = table(folder / "figure3d_spot_counts.tsv").set_index(idx + ["method"])
    compared = table(folder / "figure3d_endpoint.tsv").set_index(idx).sort_index()
    pct = (100*counts.both/counts.passing).unstack("method").sort_index()
    assert len(compared) == 384 and pct.index.equals(compared.index)
    assert np.allclose(compared[pct.columns], pct, atol=1e-12, rtol=0)
    assert (pct.source - compared.author_pct).abs().max() < 1e-12
    focused_delta = (pct.six_duplicate_groups_by_id - pct.source).abs()
    broad_delta = (pct.all_shared_ids - pct.source).abs()
    assert (focused_delta > 1e-12).sum() == 24 and focused_delta.max() < 0.0041
    assert (broad_delta > 1e-12).sum() == 356 and broad_delta.max() < 0.071
    qc = table(folder / "figure3d_qc.tsv").groupby("method").passing_spots.sum()
    assert (qc.source, qc.six_duplicate_groups_by_id, qc.all_shared_ids) == (68896, 68895, 68991)
    published = table(folder / "figure3d_published_points.tsv")
    assert len(published) == 192 and abs(published.source - published.coloc_pct).max() < 1e-10
    assert abs(published.source - published.six_duplicate_groups_by_id).gt(1e-10).sum() == 12
    tests = table(folder / "supp15_test_replay.tsv")
    decision = tests.pivot(index=["Patient_number","pair"], columns="method", values="significant")
    assert len(decision) == 192 and decision.source.equals(decision.six_duplicate_groups_by_id)
    orig = tests[tests.method.eq("source")]
    assert orig.empty_margin.sum() == 2 and (orig.loc[orig.empty_margin, "padj"] == 1).all()
    assert (orig.p_adj - orig.padj).abs().max() < 1e-10
    assert (decision.source == decision.all_shared_ids).all() and decision.source.sum() == 69
    assert abs(published.source - published.all_shared_ids).gt(1e-10).sum() == 189
    test_idx = ["Patient_number", "pair", "method"]
    comparisons = tests.set_index(test_idx)
    for key, g in counts.groupby(level=test_idx):
        contingency = np.column_stack((g.both, g.passing - g.both))
        if (contingency.sum(0) == 0).any():
            assert comparisons.at[key, "empty_margin"]
        else:
            raw_p = chi2_contingency(contingency).pvalue
            assert np.isclose(min(1, raw_p * 309 * 8), comparisons.at[key, "p_adj"], atol=1e-10, rtol=0)


def followup_replays():
    old_args, old_path = list(sys.argv), list(sys.path)
    try:
        for folder_name, program, expected in [
            ("paired_spot_control", "paired_control.py", ["per_section.tsv", "summary.json"]),
            ("targeted_pairs", "compare_pairs.py", ["pair_results.tsv", "events.tsv"]),
            ("repair_benchmark", "benchmark.py",
             ["lung_arms.tsv", "lung_changed_correct_pairs.tsv", "lung_pooled_symbols.tsv",
              "lung_wrong_pair_counterparts.tsv", "lung_slide_impact.tsv",
              "panel_arms.tsv", "validation.json"]),
        ]:
            folder = F / folder_name
            with tempfile.TemporaryDirectory() as tmp:
                script = folder / program
                sys.argv = [str(script), "--output-dir", tmp]
                sys.path.insert(0, str(folder))
                runpy.run_path(str(script), run_name="__main__")
                for name in expected:
                    assert (Path(tmp)/name).read_text(encoding="utf8") == (folder/name).read_text(encoding="utf8"), name
        folder = F / "marker_resolution"
        expected = json.loads((folder/"contrast_bounds.json").read_text())
        for folder_script in [folder/"bound_contrast.py", folder/"verify_bounds.py",
                              F/"external_order_replication/verify.py",
                              F/"independent_cohort_extension/offline_verify.py",
                              F/"prospective_holdout/verify.py",
                              F/"prospective_holdout/locus.py",
                              F/"prejoin_bounds/verify.py"]:
            sys.argv = [str(folder_script)]
            sys.path.insert(0, str(folder_script.parent))
            runpy.run_path(str(folder_script), run_name="__main__")
        assert expected == json.loads((folder/"contrast_bounds.json").read_text())
        verify_masked_donor_pairs()
    finally:
        sys.argv[:] = old_args
        sys.path[:] = old_path


def verify_masked_donor_pairs():
    blood = F / "prospective_holdout"
    official = table(blood / "official_events.tsv.gz")
    official = official.query("threshold == 1 and event == 'mismatch' and comparison != 'paired'")
    by_pair = {key: set(group["index"]) for key, group in official.groupby(["comparison", "left", "right"])}
    masked = table(F / "prejoin_bounds/pairs.tsv").fillna({"unresolved_actual_keys": ""})
    denominators = table(blood / "pair_summary.tsv").query("threshold == 1 and comparison != 'paired'")
    pair_totals = {(r.comparison, r.left, r.right): r.shared_keys for r in denominators.itertuples()}
    assert len(masked) == len(pair_totals) == 56
    for r in masked.itertuples():
        key = (r.condition, r.left, r.right)
        assert r.shared == pair_totals.pop(key)
        wrong = by_pair.pop(key, set())
        unresolved = set(r.unresolved_actual_keys.split("|")) if r.unresolved_actual_keys else set()
        assert len(wrong) == r.actual_differences and unresolved <= wrong
        assert len(wrong - unresolved) == r.certified_different and r.wrong_certifications == 0
    assert not pair_totals and not by_pair


def main():
    if not __debug__:
        raise SystemExit("Run Python without -O so all result checks are enabled.")
    sys.addaudithook(local_only)
    verify_payload()
    verify_mechanism()
    prostate()
    atlas()
    lung()
    followup_replays()

    # Render directly from the same code that generated the submitted panels.
    runpy.run_path(str(F / "v2/submission/figures_final.py"), run_name="__main__")
    fingerprints = json.loads((ROOT / "figure_hashes.json").read_text())
    for name, expected in fingerprints.items():
        path = F / "v2/submission" / name
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        assert actual == expected, (name, actual, expected)
        assert path.with_suffix(".pdf").read_bytes().startswith(b"%PDF")
    runpy.run_path(str(ROOT / "check_genejoin.py"), run_name="__main__")
    print("Source checks, study replays, three submission figures and installed genejoin passed.")


if __name__ == "__main__":
    main()
