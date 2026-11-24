#!/usr/bin/env python3
"""Extract feature tables from remote GEO ZIP files using HTTP byte ranges."""

from __future__ import annotations

import argparse
import io
import re
import socket
import urllib.request
import zipfile
from pathlib import Path


USER_AGENT = "silent-gene-identity-audit/0.1"

# The lab server's IPv6 route to NCBI is extremely slow while IPv4 is stable.
_original_getaddrinfo = socket.getaddrinfo


def _ipv4_getaddrinfo(*args, **kwargs):
    return [
        item
        for item in _original_getaddrinfo(*args, **kwargs)
        if item[0] == socket.AF_INET
    ]


socket.getaddrinfo = _ipv4_getaddrinfo


class HTTPRangeReader(io.RawIOBase):
    def __init__(self, url: str, timeout: int = 90):
        self.url = url
        self.timeout = timeout
        request = urllib.request.Request(
            url, method="HEAD", headers={"User-Agent": USER_AGENT}
        )
        with urllib.request.urlopen(request, timeout=timeout) as response:
            length = response.headers.get("Content-Length")
        if length is None:
            raise ValueError(f"No Content-Length for {url}")
        self.length = int(length)
        self.position = 0

    def readable(self) -> bool:
        return True

    def seekable(self) -> bool:
        return True

    def tell(self) -> int:
        return self.position

    def seek(self, offset: int, whence: int = io.SEEK_SET) -> int:
        if whence == io.SEEK_SET:
            position = offset
        elif whence == io.SEEK_CUR:
            position = self.position + offset
        elif whence == io.SEEK_END:
            position = self.length + offset
        else:
            raise ValueError(f"Unsupported whence: {whence}")
        if position < 0:
            raise ValueError("Negative seek position")
        self.position = position
        return position

    def read(self, size: int = -1) -> bytes:
        if self.position >= self.length:
            return b""
        if size is None or size < 0:
            end = self.length - 1
        else:
            end = min(self.length - 1, self.position + size - 1)
        request = urllib.request.Request(
            self.url,
            headers={
                "User-Agent": USER_AGENT,
                "Range": f"bytes={self.position}-{end}",
            },
        )
        with urllib.request.urlopen(request, timeout=self.timeout) as response:
            status = getattr(response, "status", response.getcode())
            content_range = response.headers.get("Content-Range")
            payload = response.read()
        if status != 206 or content_range is None:
            raise OSError(
                f"Server did not honor byte range for {self.url}: "
                f"status={status}, Content-Range={content_range!r}"
            )
        self.position += len(payload)
        return payload


def sample_bucket(sample: str) -> str:
    match = re.fullmatch(r"GSM(\d+)", sample)
    if not match:
        raise ValueError(sample)
    digits = match.group(1)
    return f"GSM{digits[:-3] or '0'}nnn"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    samples = {
        "GSM7717079": "PWH26D0",
        "GSM7717080": "PWH26D1",
        "GSM7717081": "PWH26D7",
        "GSM7717082": "PWH26D30",
        "GSM7717083": "PWH27D0",
        "GSM7717084": "PWH27D1",
        "GSM7717085": "PWH27D7",
        "GSM7717086": "PWH27D30",
        "GSM7717087": "PWH28D0",
        "GSM7717088": "PWH28D1",
        "GSM7717089": "PWH28D7",
        "GSM7717090": "PWH28D30",
    }
    for sample, label in samples.items():
        filename = f"{sample}_{label}.zip"
        url = (
            "https://ftp.ncbi.nlm.nih.gov/geo/samples/"
            f"{sample_bucket(sample)}/{sample}/suppl/{filename}"
        )
        reader = HTTPRangeReader(url)
        with zipfile.ZipFile(reader) as archive:
            candidates = [
                name
                for name in archive.namelist()
                if "__MACOSX" not in name
                and Path(name).name.lower().endswith("features.tsv.gz")
            ]
            if len(candidates) != 1:
                raise ValueError(f"{sample}: feature candidates={candidates}")
            payload = archive.read(candidates[0])
        output = args.output_dir / f"{sample}_{label}_features.tsv.gz"
        output.write_bytes(payload)
        print(sample, reader.length, candidates[0], len(payload), flush=True)


if __name__ == "__main__":
    main()
