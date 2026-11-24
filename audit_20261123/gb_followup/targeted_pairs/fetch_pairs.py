"""Download the complete human targeted/parent catalog in the pinned Scanpy loader."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import hashlib
from pathlib import Path
import re
import tempfile

import requests

HERE = Path(__file__).resolve().parent
SOURCE = "https://raw.githubusercontent.com/scverse/scanpy/1.9.8/scanpy/datasets/_datasets.py"
CACHE = Path(tempfile.gettempdir()) / "visium_target_parent_198"


def download(item):
    name, cache, expected = item
    url = f"https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/{name}/{name}_filtered_feature_bc_matrix.h5"
    path = cache / (name + ".h5")
    if not path.is_file():
        tmp = path.with_suffix(".part")
        for attempt in range(3):
            try:
                with requests.get(url, timeout=(20, 110), stream=True) as r, tmp.open("wb") as out:
                    r.raise_for_status()
                    for part in r.iter_content(1024 * 1024):
                        out.write(part)
                    if "Content-Length" in r.headers:
                        assert out.tell() == int(r.headers["Content-Length"])
                with tmp.open("rb") as f:
                    assert f.read(8) == b"\x89HDF\r\n\x1a\n"
                tmp.replace(path)
                break
            except (requests.RequestException, AssertionError):
                if attempt == 2:
                    raise
    h = hashlib.sha256(path.read_bytes()).hexdigest()
    if name in expected and h != expected[name]:
        raise ValueError(f"Published input changed: {name}")
    print(name, path.stat().st_size, flush=True)
    return dict(sample=name, url=url, bytes=path.stat().st_size, sha256=h)


def main(cache):
    r = requests.get(SOURCE, timeout=25)
    r.raise_for_status()
    loader = r.text.split("def visium_sge(", 1)[1].split("] =", 1)[0]
    names = sorted(set(re.findall("'((?:Targeted|Parent)_Visium_Human_\\w+)'", loader)))
    assert len(names) == 13 and sum(n.startswith("Targeted") for n in names) == 7
    cache.mkdir(parents=True, exist_ok=True)
    manifest = HERE / "inputs.tsv"
    expected = {}
    if manifest.is_file():
        with manifest.open(encoding="utf-8") as f:
            expected = {row["sample"]:row["sha256"] for row in csv.DictReader(f, delimiter="\t")}
    with ThreadPoolExecutor(4) as ex:
        rows = list(ex.map(download, [(name, cache, expected) for name in names]))
    with manifest.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, rows[0], delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--cache", type=Path, default=CACHE)
    main(p.parse_args().cache)
