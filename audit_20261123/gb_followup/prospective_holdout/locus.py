"""Post-lock coordinate audit of every mismatched ID pair; uses closed gene intervals."""
import argparse
import csv
import gzip
import hashlib
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from verify import HERE, table, tsv_bytes

CANONICAL = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY", "chrM"}
V19_NAME = "gencode.v19.chr_patch_hapl_scaff.annotation.gtf.gz"
V19_URL = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/" + V19_NAME


def extract(special, v19, v32):
    with gzip.open(v32, "rt", encoding="utf-8") as f:
        cache = json.load(f)
    raw = defaultdict(list)
    with gzip.open(v19, "rt", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if fields[2] != "gene":
                continue
            version = re.search(r'gene_id "([^"]+)"', fields[8])[1]
            stable = version.split(".")[0]
            if stable in special:
                raw[stable].append([fields[0], int(fields[3]), int(fields[4]), fields[6], version])
    newer = {}
    for stable in special:
        newer[stable] = [[r["chrom"], r["start"], r["end"], r["strand"], r["gene_id_version"]]
                         for r in cache["gene_records"].get(stable, [])]
    return {
        "v19": {"build": "GRCh37.p13", "url": V19_URL, "source_bytes": v19.stat().st_size,
                "source_sha256": hashlib.sha256(v19.read_bytes()).hexdigest(),
                "records": {s: raw[s] for s in sorted(special)}},
        "v32": {"build": cache["assembly"], "url": cache["source_url"],
                "source_bytes": cache["compressed_gtf_bytes"], "source_sha256": cache["compressed_gtf_sha256"],
                "records": {s: newer[s] for s in sorted(special)}}}


def classify(rows1, rows2):
    if not rows1 or not rows2:
        return "missing"
    if len(rows1) != 1 or len(rows2) != 1:
        return "non_unique"
    x, y = rows1[0], rows2[0]
    if {x[0], y[0]} - CANONICAL:
        return "noncanonical"
    if {x[0], y[0]} == {"chrX", "chrY"}:
        return "X_Y_excluded"
    if x[0] != y[0] or min(x[2], y[2]) < max(x[1], y[1]):
        return "disjoint"
    return "overlap"


def main(v19, v32):
    events = [r for r in table(HERE / "official_events.tsv.gz") if r["event"] == "mismatch"]
    feature_keys = {(r["index"], *sorted((r["left_value"], r["right_value"]))) for r in events}
    output = HERE / "locus_input.json"
    if v19 and v32:
        annotation = extract({s for _, x, y in feature_keys for s in (x, y)}, v19, v32)
        if output.exists():
            assert annotation == json.loads(output.read_text())
        else:
            output.write_text(json.dumps(annotation, indent=2) + "\n")
    else:
        annotation = json.loads(output.read_text())
    rows, aggregate, burdens = [], [], []
    sample_ids = [(r["condition"], r["sample"]) for r in table(HERE / "samples.tsv")]
    axis = table(HERE / "axis_counts.tsv.gz")
    id_rows = {r["id"]: r for r in axis}
    budgets = {(int(r["threshold"]), r["comparison"], r["sample"]): int(r["retained_UMI"])
               for r in table(HERE / "burden_by_sample.tsv") if r["event"] == "union"}
    for ref, meta in annotation.items():
        records = meta["records"]
        categories = {(key, x, y): classify(records[x], records[y]) for key, x, y in feature_keys}
        for (key, x, y), category in sorted(categories.items()):
            rows.append((ref, key, x, y, category))
        for t in (0, 1, 3, 10, "all"):
            for group in ("ctrl", "stim", "paired", "all"):
                selected = [r for r in events if (t == "all" or int(r["threshold"]) == t)
                            and (group == "all" or r["comparison"] == group)]
                strata = defaultdict(list)
                for r in selected:
                    x, y = sorted((r["left_value"], r["right_value"]))
                    strata[categories[r["index"], x, y]].append((r["index"], x, y, r["left"], r["right"]))
                if t != "all" and group != "all":
                    members = [s for c, s in sample_ids if c == group or group == "paired"]
                    hits = {(r["left"], r["left_value"]) for r in selected
                            if categories[r["index"], *sorted((r["left_value"], r["right_value"]))] == "disjoint"}
                    hits.update((r["right"], r["right_value"]) for r in selected
                                if categories[r["index"], *sorted((r["left_value"], r["right_value"]))] == "disjoint")
                    counts = sum(int(id_rows[stable][sample + "_umis"]) for sample, stable in hits)
                    burdens.append((ref, t, group, len(hits), counts, sum(budgets[t, group, s] for s in members)))
                for category in ("disjoint", "overlap", "missing", "non_unique", "noncanonical", "X_Y_excluded"):
                    group_rows = strata[category]
                    aggregate.append((ref, t, group, category, len({r[0] for r in group_rows}),
                                      len({r[1:3] for r in group_rows}), len({r[3:] for r in group_rows})))
    outputs = {
        "locus_burden.tsv": tsv_bytes(("release", "threshold", "comparison", "distinct_sample_IDs_on_disjoint_keys",
                                     "UMI_on_those_rows", "all_retained_UMI"), burdens),
        "locus_pairs.tsv": tsv_bytes(("release", "key", "id_1", "id_2", "category"), rows),
        "locus_summary.tsv": tsv_bytes(("release", "threshold", "comparison", "category", "distinct_keys",
                                       "distinct_unordered_id_pairs", "donor_sample_pairs"), aggregate)}
    for fname, data in outputs.items():
        dest = HERE / fname
        if not dest.exists() and v19 and v32:
            dest.write_bytes(data)
        assert dest.read_bytes() == data
    print("Verified", len(feature_keys), "distinct (name, unordered stable-ID pair) records in two releases")


if __name__ == "__main__":
    cmd = argparse.ArgumentParser(description=__doc__)
    cmd.add_argument("--v19", type=Path, help="official complete GENCODE v19 GTF, gzipped")
    cmd.add_argument("--v32", type=Path, help="project's GENCODE v32 primary-assembly gene cache JSON.gz")
    opt = cmd.parse_args()
    main(opt.v19, opt.v32)
