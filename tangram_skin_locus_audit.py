"""Adjudicate candidate ID mismatches by independent fixed GENCODE coordinates.

GENCODE is an external locus reference, never asserted to be the skin source GTF.
All coordinates are one-based inclusive; only nonoverlapping primary chromosome
gene intervals count toward the conservative distinct-locus endpoint.
"""
import argparse
import collections
import gzip
import hashlib
import json
from pathlib import Path
import re
import urllib.request


def gene_reference(cache, release, ids):
    cache.mkdir(parents=True, exist_ok=True)
    name = f"gencode.v{release}.primary_assembly.annotation.gtf.gz"
    url = f"https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{release}/{name}"
    meta_path = cache / f"gencode.v{release}.gene_metadata.json.gz"
    if meta_path.exists():
        with gzip.open(meta_path, "rt") as f:
            meta = json.load(f)
    else:
        raw_path = cache / name
        if not raw_path.exists():
            with urllib.request.urlopen(url, timeout=120) as r, raw_path.open("wb") as f:
                while chunk := r.read(1024 * 1024):
                    f.write(chunk)
        sha = hashlib.sha256()
        with raw_path.open("rb") as f:
            while chunk := f.read(1024 * 1024):
                sha.update(chunk)
        genes = collections.defaultdict(list)
        headers = []
        with gzip.open(raw_path, "rt") as f:
            for line in f:
                if line.startswith("#"):
                    headers.append(line.rstrip())
                    continue
                cols = line.rstrip("\n").split("\t")
                if len(cols) != 9 or cols[2] != "gene":
                    continue
                attrs = dict(re.findall(r'(\w+) "([^"]+)"', cols[8]))
                full = attrs["gene_id"]
                stable = full.split(".")[0]
                genes[stable].append({"gene_id": stable, "gene_id_version": full,
                                      "gene_name": attrs.get("gene_name"),
                                      "gene_type": attrs.get("gene_type"),
                                      "chrom": cols[0], "start": int(cols[3]),
                                      "end": int(cols[4]), "strand": cols[6]})
        meta = {"source_url": url, "release": release, "assembly": "GRCh38",
                "compressed_gtf_bytes": raw_path.stat().st_size, "compressed_gtf_sha256": sha.hexdigest(),
                "gtf_headers": headers, "gene_records": dict(genes)}
        with gzip.open(meta_path, "wt") as f:
            json.dump(meta, f, separators=(",", ":"))
    subset = {k: meta["gene_records"].get(k, []) for k in sorted(ids)}
    return {k: v for k, v in meta.items() if k != "gene_records"} | {
        "gene_metadata_cache": meta_path.name, "n_genes_in_cache": len(meta["gene_records"]),
        "candidate_records": subset, "candidate_ids_missing": [k for k, v in subset.items() if not v]}


def locus_class(a, b):
    if not a or not b:
        return {"classification": "unresolved_missing_id", "strict_distinct_locus": False}
    canonical = re.compile(r"^chr(?:[1-9]|1[0-9]|2[0-2]|X|Y|M)$")
    if len(a) != 1 or len(b) != 1 or not all(canonical.match(x["chrom"]) for x in a + b):
        return {"classification": "unresolved_noncanonical_or_multilocus", "strict_distinct_locus": False}
    x, y = a[0], b[0]
    if x["chrom"] != y["chrom"]:
        # Conservatively exclude potential PAR X/Y equivalents.
        if {x["chrom"], y["chrom"]} == {"chrX", "chrY"}:
            return {"classification": "unresolved_xy_equivalence", "strict_distinct_locus": False}
        return {"classification": "distinct_chromosomes", "strict_distinct_locus": True, "gap_bp": None}
    overlap = max(0, min(x["end"], y["end"]) - max(x["start"], y["start"]) + 1)
    if overlap:
        return {"classification": "overlapping_gene_intervals_excluded", "strict_distinct_locus": False,
                "overlap_bp": overlap, "same_strand": x["strand"] == y["strand"],
                "interpretation": "Conservatively excluded from distinct-locus endpoint; overlap alone does not prove biological identity."}
    gap = max(x["start"], y["start"]) - min(x["end"], y["end"]) - 1
    return {"classification": "nonoverlapping_same_chromosome", "strict_distinct_locus": True,
            "gap_bp": gap, "same_strand": x["strand"] == y["strand"]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", type=Path, required=True)
    ap.add_argument("--cache", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--release", type=int, default=32)
    args = ap.parse_args()
    raw = json.loads(args.input.read_text(encoding="utf8"))
    candidates = collections.defaultdict(set)
    for p in raw["pairs"]:
        for c in p["conflicts"]:
            candidates[c["key"]].add(tuple(sorted([c["reference_id"], c["spatial_id"]])))
    ids = {i for pairs in candidates.values() for pair in pairs for i in pair}
    ref = gene_reference(args.cache, args.release, ids)
    adjudications = []
    lookup = {}
    for key, pairs in sorted(candidates.items()):
        for a, b in sorted(pairs):
            x, y = ref["candidate_records"][a], ref["candidate_records"][b]
            row = {"key": key, "id_a": a, "id_b": b, "locus_a": x, "locus_b": y, **locus_class(x, y)}
            adjudications.append(row)
            lookup[tuple(sorted([a, b]))] = row
    pairs_out = []
    groups = collections.defaultdict(list)
    pooled_reference_features = {}
    for p in raw["pairs"]:
        strict = [c for c in p["conflicts"] if lookup[tuple(sorted([c["reference_id"], c["spatial_id"]]))]["strict_distinct_locus"]]
        training = [c for c in strict if c["enters_default_training"]]
        row = {"reference": p["reference"], "spatial": p["spatial"], "mode": p["mode"],
               "strict_overlap_conflicts": len(strict), "strict_training_conflicts": len(training),
               "strict_training_keys": [c["key"] for c in training],
               "strict_training_reference_counts": sum(c["reference_count"] for c in training),
               "strict_training_spatial_counts": sum(c["spatial_count"] for c in training),
               "reference_total_counts": p["reference_total_counts"],
               "spatial_total_counts": p["spatial_total_counts"]}
        row["strict_reference_fraction"] = row["strict_training_reference_counts"] / row["reference_total_counts"]
        row["strict_spatial_fraction"] = row["strict_training_spatial_counts"] / row["spatial_total_counts"]
        pairs_out.append(row)
        scope = "pooled" if str(p["reference"]).startswith("POOLED") or "pooled" in str(p["reference"]).lower() else "individual"
        groups[(p["mode"], scope)].append(row)
        if p["mode"] == "filter_then_unique" and scope == "pooled":
            for c in training:
                if c["reference_id"] in pooled_reference_features:
                    assert pooled_reference_features[c["reference_id"]] == c["reference_count"]
                pooled_reference_features[c["reference_id"]] = c["reference_count"]
    summary = [{"mode": mode, "reference_scope": scope, "sample_pairs": len(ps),
                "strict_positive_training_pairs": sum(p["strict_training_conflicts"] > 0 for p in ps),
                "max_strict_training_conflicts": max(p["strict_training_conflicts"] for p in ps),
                "strict_training_keys": sorted({k for p in ps for k in p["strict_training_keys"]})}
               for (mode, scope), ps in groups.items()]
    ps = groups.get(("filter_then_unique", "pooled"), [])
    assert len({p["spatial"] for p in ps}) == len(ps), "Spatial UMI aggregate must not double-count samples"
    pos = [p for p in ps if p["strict_training_conflicts"]]
    burden = {"scope": "Pooled reference versus each spatial sample exactly once; no sum across the dependent 192 individual-reference pairs.",
              "spatial_sample_count": len(ps), "spatial_total_umis": sum(p["spatial_total_counts"] for p in ps),
              "strict_affected_spatial_umis": sum(p["strict_training_spatial_counts"] for p in ps),
              "pooled_reference_total_umis": ps[0]["reference_total_counts"] if ps else 0,
              "pooled_reference_unique_affected_feature_umis": sum(pooled_reference_features.values()),
              "pooled_reference_unique_affected_features": pooled_reference_features}
    burden["strict_affected_spatial_umi_fraction"] = burden["strict_affected_spatial_umis"] / burden["spatial_total_umis"]
    for label, subset in [("all16", ps), ("positive11", pos)]:
        burden[label] = {field: [min(p[field] for p in subset), max(p[field] for p in subset)]
                         for field in ["strict_training_reference_counts", "strict_training_spatial_counts", "strict_reference_fraction", "strict_spatial_fraction"]}
    result = {"analysis_date": "2026-09-19", "input_file": args.input.name,
              "input_sha256": hashlib.sha256(args.input.read_bytes()).hexdigest(),
              "scope": "Independent coordinate adjudication of counterfactual Tangram feature identity conflicts, not original paper error adjudication.",
              "reference_provenance_limit": "Official GENCODE fixed release is an external locus reference. We have not established that this is the exact GTF used to build the deposited skin feature table.",
              "strict_rule": "Both IDs resolve uniquely to nonoverlapping canonical-chromosome gene intervals in the same fixed GENCODE release. All overlapping or missing/multilocus/XY-equivalence cases excluded.",
              "reference": ref, "candidate_id_pairs": adjudications, "summary": summary,
              "pooled_burden_without_spatial_double_counting": burden, "pairs": pairs_out}
    args.output.write_text(json.dumps(result, indent=2), encoding="utf8")
    print(json.dumps({"missing_ids": ref["candidate_ids_missing"], "candidates": [{"key": r["key"], "classification": r["classification"]} for r in adjudications], "summary": summary}, indent=2), flush=True)


if __name__ == "__main__":
    main()
