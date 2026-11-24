"""Independently count default-training identity conflicts and check tool-specific controls.

Run with Python 3.9+; pass --write to refresh order_controls.json. Inputs are
the archived per-pair and per-feature measurements and fixed GENCODE v32 loci.
The original raw-count reconstruction is audit_tangram_skin_real_counts.py.
"""
import argparse
from collections import defaultdict
import gzip
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
CANONICAL = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY", "chrM"}


def read(name):
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def distinct_loci(genes, left, right):
    a, b = genes.get(left, []), genes.get(right, [])
    if len(a) != 1 or len(b) != 1:
        return False
    x, y = a[0], b[0]
    chromosomes = {x["chrom"], y["chrom"]}
    if not chromosomes <= CANONICAL or chromosomes == {"chrX", "chrY"}:
        return False
    return len(chromosomes) == 2 or min(x["end"], y["end"]) < max(x["start"], y["start"])


def skin():
    source = read("tangram_skin_real_counts_result.json")
    with gzip.open(ROOT / "gencode.v32.gene_metadata.json.gz", "rt", encoding="utf-8") as stream:
        annotation = json.load(stream)
    assert annotation["release"] == 32
    sources = source["sources"]
    assert len(sources) == len({s["sample"] for s in sources}) == 28
    assert {s["n_genes"] for s in sources} == {36600}
    groups = defaultdict(list)
    nc_umis = []
    keys = set()
    for row in source["pairs"]:
        scope = "pooled" if row["reference"].startswith("POOLED_") else "individual"
        hits = [c for c in row["conflicts"] if c["enters_default_training"]]
        assert all(c["reference_id"] != c["spatial_id"] and
                   min(c["reference_count"], c["spatial_count"]) > 0 for c in hits)
        strict = [c for c in hits if distinct_loci(annotation["gene_records"], c["reference_id"], c["spatial_id"])]
        groups[row["mode"], scope].append((row["reference"], row["spatial"], len(hits), len(strict)))
        if row["mode"] == "filter_then_unique" and scope == "pooled":
            nc_umis.append((row["spatial"], sum(c["spatial_count"] for c in strict), row["spatial_total_counts"]))
        keys.update(c["key"] for c in strict)
    comparisons = []
    for (mode, scope), rows in sorted(groups.items()):
        assert len(rows) == len({(r[0], r[1]) for r in rows}) == (16 if scope == "pooled" else 192)
        comparisons.append(dict(order=mode, reference=scope, pairs=len(rows),
                                different_id_pairs=sum(r[2] > 0 for r in rows),
                                disjoint_locus_pairs=sum(r[3] > 0 for r in rows),
                                maximum_disjoint_keys_in_one_pair=max(r[3] for r in rows)))
    assert len(comparisons) == 6 and len(nc_umis) == len({x[0] for x in nc_umis}) == 16
    checks = source["source_function_checks"]
    src_sha = hashlib.sha256((ROOT / "tangram_mapping_utils_pinned.py").read_bytes()).hexdigest()
    for c in checks:
        scope = "pooled" if c["reference"].startswith("POOLED_") else "individual"
        hit, = [r for r in groups[c["mode"], scope] if r[:2] == (c["reference"], c["spatial"])]
        assert c["training_conflicts"] == hit[2] and c["tangram_source_sha256"] == src_sha
        assert c["naming_filter_and_training_set_assertions"] == "passed"
    assert len(checks) == 9 and len({c["tangram_commit"] for c in checks}) == 1
    return dict(scope="Counterfactual Tangram preprocessing on deposited skin counts; no optimization or published skin endpoint was tested.",
                ordered_feature_table_sha256=source["raw_ordered_feature_sha256"], gene_rows=36600,
                comparisons=comparisons, disjoint_training_keys=sorted(keys),
                pooled_distinct_spatial_umi=sum(x[1] for x in nc_umis),
                pooled_spatial_total_umi=sum(x[2] for x in nc_umis),
                gencode_v32_source_sha256=annotation["compressed_gtf_sha256"],
                pinned_pp_adatas_commit=checks[0]["tangram_commit"], source_replay_checks=len(checks))


def scvelo():
    count = read("scvelo_10x_count_burden.json")
    replay = read("scvelo_real_10x_merge_reproduction.json")
    rows = count["affected_genes"]
    assert all(r["rna_id"] != r["loom_id"] for r in rows)
    for row in rows:
        loci = row["loci"]
        assert len(set(loci["Chromosome"])) == 1
        assert min(map(int, loci["End"])) < max(map(int, loci["Start"]))
        in_selection = row["tutorial_neutrophils"]
        assert in_selection["rna_count"] == 0
        assert all(in_selection[layer]["wrong_count"] == in_selection[layer]["correct_count"] == 0
                   for layer in ("spliced", "unspliced", "ambiguous"))
    assert replay["full_input_ID_to_symbol_maps_identical"]
    selected = [r for r in replay["results"] if r["cohort"] == "tutorial_neutrophils"]
    assert len(selected) == 2
    for r in selected:
        assert r["cells"] == count["tutorial_selected_cells"]
        assert r["identity_rebindings"] == (len(rows) if r["matching_key"] == "locally_unique_symbol" else 0)
    return dict(scope="Observed input row-order difference in a published 10x tutorial; all affected rows zero in selected cells.",
                source_merge_commit=replay["source_commit"], primary_tutorial_selected_cells=count["tutorial_selected_cells"],
                primary_tutorial_affected_umi=0, switched_keys=sorted(r["name"] for r in rows))


def signac():
    source = read("signac_real_identity_audit.json")
    hvg = [r for r in source["conflicts"] if r["in_published_rna_hvg"]]
    duplicate = [r for r in hvg if max(r["rna_symbol_group_size"], r["atac_symbol_group_size"]) > 1]
    assert len(duplicate) == 1 and duplicate[0]["key"] == "ECSCR" and duplicate[0]["same_or_overlapping_locus"]
    return dict(scope="Gene-activity annotation candidates against the published RNA HVGs; fragments and anchor output untested.",
                published_rna_hvg=source["dimensions"]["published_rna_hvg"],
                mismatched_candidate_hvg_keys=sorted(r["key"] for r in hvg),
                duplicate_associated_keys=[r["key"] for r in duplicate],
                disjoint_duplicate_associated_keys=[r["key"] for r in duplicate if not r["same_or_overlapping_locus"]])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true", help="replace output with the recalculated evidence")
    args = parser.parse_args()
    result = {"tangram": skin(), "scvelo": scvelo(), "signac": signac()}
    path = HERE / "order_controls.json"
    if args.write:
        path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    else:
        assert json.loads(path.read_text(encoding="utf-8")) == result, "stored results differ"
    print("Verified 624 skin comparisons, 9 source-function checks and the scVelo and Signac controls.")


if __name__ == "__main__":
    main()
