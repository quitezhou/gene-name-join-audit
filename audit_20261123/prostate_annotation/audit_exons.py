"""Measure genomic exon overlap of the ten prostate feature-ID pairs (GENCODE v32).

Usage: python audit_exons.py PATH_TO_gencode.v32.primary_assembly.annotation.gtf.gz
"""
import argparse
import csv
import gzip
import hashlib
import json
import re
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
SOURCE = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.primary_assembly.annotation.gtf.gz"
SHA = "0c36c37abb0c88637b04a8e401c93f0199a65bbdd0702c9e13bdfbadc053fe1e"


def exon_union(intervals):
    """Merge closed intervals across all annotated transcripts of one gene."""
    merged = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(end, merged[-1][1]))
        else:
            merged.append((start, end))
    return merged


def audit(path):
    if hashlib.sha256(path.read_bytes()).hexdigest() != SHA:
        raise ValueError("GTF checksum differs from frozen GENCODE release")
    source = HERE.parents[1] / "prostate_spatial_versions_result.json"
    pairs = json.loads(source.read_text(encoding="utf-8"))["comparison_after_per_table_anndata_uniquification"]["mismatch_examples"]
    ids = {pair[k] for pair in pairs for k in ["first_stable_id", "second_stable_id"]}
    exon = defaultdict(list)
    genes = {}
    with gzip.open(path, "rt", encoding="utf-8") as data:
        for line in data:
            if line.startswith("#"):
                continue
            chrom, _, kind, start, end, _, strand, _, attrs = line.rstrip().split("\t")
            if kind not in ("gene", "exon"):
                continue
            full = re.search(r'gene_id "([^"]+)"', attrs)[1]
            gene_id = full.split(".")[0]
            if gene_id not in ids:
                continue
            if kind == "exon":
                exon[gene_id].append((int(start), int(end)))
            else:
                if gene_id in genes:
                    raise ValueError(f"Multiple coordinates for {gene_id}")
                biotype = re.search(r'gene_type "([^"]+)"', attrs)[1]
                genes[gene_id] = (chrom, strand, full, biotype)
    rows = []
    for pair in pairs:
        first, second = pair["first_stable_id"], pair["second_stable_id"]
        row = dict(label=pair["local_unique_name"], old_table_id=first, new_table_id=second)
        for key, gid in [("old", first), ("new", second)]:
            chrom, strand, version, biotype = genes.get(gid, ("", "", "", ""))
            row.update({key + "_chrom": chrom, key + "_strand": strand,
                        key + "_version": version, key + "_biotype": biotype})
        if first not in genes or second not in genes:
            row["class_v32"] = "unresolved_missing_annotation"
        else:
            left, right = exon_union(exon[first]), exon_union(exon[second])
            if not left or not right:
                raise ValueError(f"Gene without exons: {row['label']}")
            same_chr = genes[first][0] == genes[second][0]
            common = sum(max(0, min(b, d) - max(a, c) + 1)
                         for a, b in left for c, d in right) if same_chr else 0
            row.update(old_exonic_bp=sum(b - a + 1 for a, b in left),
                       new_exonic_bp=sum(b - a + 1 for a, b in right), shared_exonic_bp=common)
            if common:
                row["class_v32"] = ("shared_exons_same_strand" if genes[first][1] == genes[second][1]
                                    else "shared_exons_opposite_strand")
            else:
                row["class_v32"] = "no_shared_exonic_bases"
        rows.append(row)
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with (HERE / "exon_overlap.tsv").open("w", newline="", encoding="utf-8") as output:
        writer = csv.DictWriter(output, fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return rows


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gtf", type=Path, help=SOURCE)
    audit(parser.parse_args().gtf)
