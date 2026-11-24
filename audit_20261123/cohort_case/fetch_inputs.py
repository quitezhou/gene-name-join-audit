"""Fetch verifiable published feature tables and inventory paths. Run with Python and requests."""
import csv
import hashlib
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import requests

OUT = Path(__file__).resolve().parent
SOURCE = OUT / "source"
SOURCE.mkdir(exist_ok=True)


def get(url):
    response = requests.get(url, timeout=(15, 90))
    response.raise_for_status()
    return response


def fetch_file(entry):
    out = SOURCE / (entry["file_id"] + ".tsv.gz")
    if not out.exists():
        data = get(entry["url"]).content
        if hashlib.sha256(data).hexdigest() != entry["sha256"]:
            raise ValueError(f"checksum mismatch: {entry['path']}")
        out.write_bytes(data)
    if hashlib.sha256(out.read_bytes()).hexdigest() != entry["sha256"]:
        raise ValueError(f"local checksum mismatch: {entry['path']}")


def write(name, entries):
    with (OUT / name).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, list(entries[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(entries)


def main():
    inventory = []
    for dataset, version in [("y3pxvr4xbf", 2), ("3bs5f8mvbs", 1)]:
        url = f"https://data.mendeley.com/public-api/datasets/{dataset}"
        record = get(url).json()
        if int(record["version"]) != version:
            raise ValueError(f"Unexpected dataset version for {dataset}")
        folders = {f["id"]: f for f in get(f"{url}/folders/{version}").json()}
        for f in record["files"]:
            components = [f["filename"]]
            folder = f.get("folder_id")
            while folder:
                components.insert(0, folders[folder]["name"])
                folder = folders[folder].get("parent_id")
            detail = f["content_details"]
            inventory.append(dict(dataset=dataset, version=version, path="/".join(components),
                                  file_id=f["id"], size=detail["size"], sha256=detail["sha256_hash"],
                                  url=detail["download_url"]))
    inventory.sort(key=lambda row: (row["dataset"],row["path"]))
    write("deposit_inventory.tsv", inventory)
    selected = [r for r in inventory if r["path"].endswith("features.tsv.gz")]
    with ThreadPoolExecutor(max_workers=6) as pool:
        list(pool.map(fetch_file, selected))
    print(f"{len(selected)} verified files across {len(inventory)} deposit entries")

if __name__ == "__main__":
    main()
