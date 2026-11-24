"""Retrieve a deposit entry by exact relative path into the local temporary cache."""
import argparse
import csv
import hashlib
import tempfile
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor

import requests

OUT = Path(__file__).resolve().parent


def cache(r):
    path = Path(tempfile.gettempdir()) / ("hwat_" + r["file_id"] + "_" + Path(r["path"]).name)
    if not path.exists():
        tmp = path.with_suffix(path.suffix + ".part")
        with requests.get(r["url"], stream=True, timeout=(20, 150)) as resp:
            resp.raise_for_status()
            with tmp.open("wb") as file:
                for block in resp.iter_content(1024*1024):
                    file.write(block)
        with tmp.open("rb") as file:
            if hashlib.file_digest(file, "sha256").hexdigest() != r["sha256"]:
                raise ValueError(r["path"])
        tmp.replace(path)
    with path.open("rb") as file:
        if hashlib.file_digest(file, "sha256").hexdigest() != r["sha256"]:
            raise ValueError(r["path"])
    return path


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("paths", nargs="+", help="exact paths as listed in deposit_inventory.tsv")
    args = parser.parse_args()
    with (OUT/"deposit_inventory.tsv").open(encoding="utf8") as src:
        entries = {r["path"]:r for r in csv.DictReader(src, delimiter="\t")}
    with ThreadPoolExecutor(max_workers=3) as pool:
        for path in pool.map(cache, [entries[s] for s in args.paths]):
            print(path, flush=True)
