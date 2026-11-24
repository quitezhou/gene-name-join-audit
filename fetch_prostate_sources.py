"""Retrieve checksummed publisher/author sources for the optional prostate refit."""
import argparse
import hashlib
import json
from pathlib import Path
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parent
CASE = ROOT / "audit_20261123/sentinel_case"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--all", action="store_true", help="Also retrieve the original notebooks and other workbooks")
    args = parser.parse_args()
    folder = CASE / "source"
    folder.mkdir(exist_ok=True)
    for name, rec in json.loads((CASE / "primary_sources.json").read_text()).items():
        if not args.all and name != "MOESM13.xlsx":
            continue
        path = folder / name
        if path.exists():
            data = path.read_bytes()
        else:
            with urlopen(rec["url"], timeout=90) as src:
                data = src.read()
        if hashlib.sha256(data).hexdigest() != rec["sha256"] or len(data) != rec["bytes"]:
            raise ValueError(f"Source mismatch for {name}")
        if not path.exists():
            path.write_bytes(data)
        print(path)


if __name__ == "__main__":
    main()
