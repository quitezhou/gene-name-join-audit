"""Retrieve public GSE241132 count members; summarize with the stdlib only.

Each member is downloaded using the ZIP directory's HTTP byte ranges and CRC
checked by zipfile. No sequencing reads or controlled-access files are used.
"""
from __future__ import annotations
import argparse
import gzip
import hashlib
import io
import json
from pathlib import Path
import zipfile

from extract_remote_zip_features import HTTPRangeReader, sample_bucket

SAMPLES = {
    "GSM7717079": "PWH26D0", "GSM7717080": "PWH26D1",
    "GSM7717081": "PWH26D7", "GSM7717082": "PWH26D30",
    "GSM7717083": "PWH27D0", "GSM7717084": "PWH27D1",
    "GSM7717085": "PWH27D7", "GSM7717086": "PWH27D30",
    "GSM7717087": "PWH28D0", "GSM7717088": "PWH28D1",
    "GSM7717089": "PWH28D7", "GSM7717090": "PWH28D30",
}


def get_member(zf, suffix, cache):
    found = [x for x in zf.infolist() if not x.is_dir() and x.filename.endswith(suffix)
             and not x.filename.startswith("__MACOSX/")]
    if len(found) != 1:
        raise ValueError(f"Expected one {suffix}, got {[x.filename for x in found]}")
    info = found[0]
    path = cache / suffix
    if path.exists():
        blob = path.read_bytes()
        import zlib
        if len(blob) != info.file_size or zlib.crc32(blob) != info.CRC:
            raise ValueError(f"Cached ZIP member integrity mismatch: {path}")
    else:
        print("FETCH", info.filename, info.compress_size, flush=True)
        blob = zf.read(info)
        path.write_bytes(blob)
    return blob, {"member": info.filename, "size": len(blob), "crc32": info.CRC,
                  "sha256": hashlib.sha256(blob).hexdigest()}


def extract(sample, label, root):
    target = root / f"{sample}_count_summary.json.gz"
    if target.exists():
        with gzip.open(target, "rt") as stream:
            prior = json.load(stream)
        print("CACHED", sample, prior["nnz"], flush=True)
        return
    cache = root / "members" / sample
    cache.mkdir(parents=True, exist_ok=True)
    url = (f"https://ftp.ncbi.nlm.nih.gov/geo/samples/{sample_bucket(sample)}/"
           f"{sample}/suppl/{sample}_{label}.zip")
    reader = HTTPRangeReader(url, timeout=180)
    with zipfile.ZipFile(reader) as zf:
        print("ARCHIVE", sample, reader.length, [x.filename for x in zf.infolist()], flush=True)
        features_blob, fmeta = get_member(zf, "features.tsv.gz", cache)
        matrix_blob, mmeta = get_member(zf, "matrix.mtx.gz", cache)
        barcode_blob, bmeta = get_member(zf, "barcodes.tsv.gz", cache)
    features = [line.split("\t") for line in gzip.decompress(features_blob).decode().splitlines()]
    barcode_count = len(gzip.decompress(barcode_blob).splitlines())
    counts = [0] * len(features)
    positive_entries = [0] * len(features)
    observed = 0
    with gzip.GzipFile(fileobj=io.BytesIO(matrix_blob)) as stream:
        header = stream.readline().decode().strip()
        if header != "%%MatrixMarket matrix coordinate integer general":
            raise ValueError(f"Unexpected MatrixMarket schema: {header}")
        line = stream.readline()
        while line.startswith(b"%"):
            line = stream.readline()
        n_genes, n_cells, nnz = map(int, line.split())
        assert n_genes == len(features) and n_cells == barcode_count
        for line in stream:
            if not line.strip():
                continue
            gene, cell, value = map(int, line.split())
            assert 1 <= gene <= n_genes and 1 <= cell <= n_cells and value >= 0
            counts[gene - 1] += value
            positive_entries[gene - 1] += value > 0
            observed += 1
        assert observed == nnz
    result = {
        "series": "GSE241132", "sample": sample, "label": label,
        "source_url": url, "archive_bytes": reader.length,
        "features": fmeta, "matrix": mmeta, "barcodes": bmeta,
        "n_genes": n_genes, "n_cells": n_cells, "nnz": nnz,
        "total_counts": sum(counts),
        "feature_rows": [dict(stable_id=f[0], symbol=f[1],
                              feature_type=f[2] if len(f)>2 else "Gene Expression",
                              total_count=c, positive_entries=n)
                         for f,c,n in zip(features,counts,positive_entries)],
        "interpretation": "Per-feature sums of actual deposited nonnegative counts. positive_entries are positive coordinate entries, not a deduplicated cell count; >0 exactly identifies genes surviving min_cells=1.",
    }
    with gzip.open(target, "wt", encoding="utf-8") as stream:
        json.dump(result, stream, ensure_ascii=False)
    print("COMPLETE", sample, n_genes, n_cells, nnz, sum(counts), flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--samples", nargs="*", default=list(SAMPLES))
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for sample in args.samples:
        extract(sample, SAMPLES[sample], args.output_dir)


if __name__ == "__main__":
    main()
