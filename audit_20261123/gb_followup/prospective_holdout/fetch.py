"""Fetch the six unmodified GEO deposit files and record their byte-level SHA-256."""
import argparse
import concurrent.futures
import hashlib
import json
from pathlib import Path
import requests

SERIES = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE96nnn/GSE96583/suppl/"
SAMPLE = "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2560nnn/{}/suppl/"
URLS = [SERIES + "GSE96583_batch2." + suffix for suffix in (
    "genes.tsv.gz", "total.tsne.df.tsv.gz")]
for gsm, matrix in (("GSM2560248", "2.1"), ("GSM2560249", "2.2")):
    for suffix in ("barcodes.tsv.gz", matrix + ".mtx.gz"):
        URLS.append(SAMPLE.format(gsm) + gsm + "_" + suffix)


def download(url, target):
    path = target / url.rsplit("/", 1)[1]
    if not path.exists():
        temporary = path.with_suffix(".part")
        with requests.get(url, timeout=(15, 180), stream=True) as response:
            response.raise_for_status()
            with temporary.open("wb") as out:
                for block in response.iter_content(512 * 1024):
                    out.write(block)
        temporary.replace(path)
    with path.open("rb") as fin:
        digest = hashlib.file_digest(fin, "sha256").hexdigest()
    return {"name": path.name, "url": url, "bytes": path.stat().st_size, "sha256": digest}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cache", type=Path)
    cli = parser.parse_args()
    cli.cache.mkdir(exist_ok=True, parents=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=6) as pool:
        files = list(pool.map(lambda url: download(url, cli.cache), URLS))
    expected = Path(__file__).with_name("sources.json")
    if expected.exists() and files != json.loads(expected.read_text()):
        raise ValueError("Input hash mismatch")
    expected.write_text(json.dumps(files, indent=2) + "\n")
