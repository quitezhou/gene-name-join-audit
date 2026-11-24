"""Verify or fetch the pinned study files, optionally regenerate RNA exports."""
import argparse
import csv
import hashlib
from pathlib import Path
import subprocess
import sys

import requests

HERE = Path(__file__).resolve().parent
CASE = HERE.parents[1] / "cohort_case"
sys.path.insert(0, str(CASE))
from cache_entry import cache


def verify(file, expected):
    with file.open("rb") as src:
        observed = hashlib.file_digest(src, "sha256").hexdigest()
    if observed != expected:
        raise ValueError(f"Digest mismatch: {file}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rscript", help="Also re-export four deposited RDS objects with this Rscript")
    args = parser.parse_args()
    target = HERE / "source"
    target.mkdir(exist_ok=True)
    with (HERE / "source_manifest.tsv").open() as src:
        for row in csv.DictReader(src, delimiter="\t"):
            file = target / row["file"]
            if not file.exists():
                staging = file.with_suffix(file.suffix + ".part")
                with requests.get(row["url"], stream=True, timeout=(20, 120)) as response:
                    response.raise_for_status()
                    with staging.open("wb") as out:
                        for chunk in response.iter_content(1 << 20):
                            out.write(chunk)
                verify(staging, row["sha256"])
                staging.replace(file)
            else:
                verify(file, row["sha256"])
    if args.rscript:
        with (CASE / "deposit_inventory.tsv").open() as src:
            inv = {r["path"]: r for r in csv.DictReader(src, delimiter="\t")}
        for name in ["ly", "B", "endo", "scFAP"]:
            key = "Source_Data_complete/Source Data II_Seurat Files/" + name + ".rds"
            subprocess.run([args.rscript, str(HERE / "export_marker_vectors.R"),
                            str(cache(inv[key])), str(target / (name + "_rows.tsv.gz"))], check=True)


if __name__ == "__main__":
    main()
