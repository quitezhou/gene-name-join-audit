"""Download and fingerprint every GSE278936 feature table from public GEO."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
import time
import urllib.request
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


def parse_filelist(path: Path) -> list[tuple[str, int]]:
    result = []
    for line in path.read_text(encoding="utf-8").splitlines():
        fields = line.split("\t")
        if len(fields) >= 5 and fields[1].endswith("_features.tsv.gz"):
            result.append((fields[1], int(fields[3])))
    return result


def geo_url(filename: str) -> str:
    match = re.match(r"(GSM\d+)_", filename)
    if not match:
        raise ValueError(f"No GSM accession in {filename}")
    accession = match.group(1)
    bucket = accession[:-3] + "nnn"
    return (
        f"https://ftp.ncbi.nlm.nih.gov/geo/samples/{bucket}/{accession}/suppl/"
        f"{filename}"
    )


def download_one(filename: str, expected_bytes: int, output_dir: Path) -> dict:
    destination = output_dir / filename
    url = geo_url(filename)
    error = None
    for attempt in range(1, 13):
        try:
            if not destination.exists() or destination.stat().st_size != expected_bytes:
                with urllib.request.urlopen(url, timeout=20) as response:
                    payload = response.read()
                if len(payload) != expected_bytes:
                    raise ValueError(
                        f"expected {expected_bytes} bytes, received {len(payload)}"
                    )
                destination.write_bytes(payload)
            compressed = destination.read_bytes()
            decompressed = gzip.decompress(compressed)
            rows = sum(bool(line) for line in decompressed.splitlines())
            return {
                "file": filename,
                "expected_compressed_bytes": expected_bytes,
                "compressed_sha256": hashlib.sha256(compressed).hexdigest(),
                "decompressed_sha256": hashlib.sha256(decompressed).hexdigest(),
                "rows": rows,
                "url": url,
            }
        except Exception as exc:  # network retries are reported in the final JSON
            error = f"{type(exc).__name__}: {exc}"
            if destination.exists() and destination.stat().st_size != expected_bytes:
                destination.unlink()
            time.sleep(min(2 * attempt, 10))
    raise RuntimeError(f"Failed {filename}: {error}")


def fingerprint_existing(filename: str, expected_bytes: int, path: Path) -> dict:
    compressed = path.read_bytes()
    if len(compressed) != expected_bytes:
        raise ValueError(
            f"{filename}: expected {expected_bytes} bytes, found {len(compressed)}"
        )
    decompressed = gzip.decompress(compressed)
    return {
        "file": filename,
        "expected_compressed_bytes": expected_bytes,
        "compressed_sha256": hashlib.sha256(compressed).hexdigest(),
        "decompressed_sha256": hashlib.sha256(decompressed).hexdigest(),
        "rows": sum(bool(line) for line in decompressed.splitlines()),
        "url": geo_url(filename),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--filelist", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=12)
    parser.add_argument("--offline", action="store_true")
    parser.add_argument("--extra-dir", type=Path)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    listed = parse_filelist(args.filelist)
    completed = []
    failures = []
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {}
        for filename, size in listed:
            if args.offline:
                candidates = [args.output_dir / filename]
                if args.extra_dir is not None:
                    candidates.append(args.extra_dir / filename)
                existing = next((path for path in candidates if path.exists()), None)
                if existing is None:
                    failures.append({"file": filename, "error": "missing in offline mode"})
                    continue
                future = executor.submit(fingerprint_existing, filename, size, existing)
            else:
                future = executor.submit(download_one, filename, size, args.output_dir)
            futures[future] = filename
        for future in as_completed(futures):
            filename = futures[future]
            try:
                completed.append(future.result())
                print(f"complete {filename}", flush=True)
            except Exception as exc:
                failures.append({"file": filename, "error": str(exc)})
                print(f"failed {filename}: {exc}", flush=True)

    completed.sort(key=lambda row: row["file"])
    fingerprint_counts = Counter(row["decompressed_sha256"] for row in completed)
    cohorts = []
    for fingerprint, count in sorted(fingerprint_counts.items()):
        members = [row for row in completed if row["decompressed_sha256"] == fingerprint]
        cohorts.append(
            {
                "decompressed_sha256": fingerprint,
                "files": count,
                "rows": sorted({row["rows"] for row in members}),
                "compressed_byte_sizes": sorted(
                    {row["expected_compressed_bytes"] for row in members}
                ),
                "first_file": members[0]["file"],
                "last_file": members[-1]["file"],
            }
        )
    result = {
        "analysis_date": "2026-09-19",
        "listed_feature_files": len(listed),
        "completed_feature_files": len(completed),
        "failures": failures,
        "exact_decompressed_content_cohorts": cohorts,
        "files": completed,
    }
    args.output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps({key: result[key] for key in result if key != "files"}, indent=2))
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
