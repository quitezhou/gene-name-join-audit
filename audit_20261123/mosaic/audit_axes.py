"""Audit released MOSAIC Chromium axes against a public 10x 2020-A feature table.

Run with Python, h5py, fsspec, aiohttp, requests. No expression or donor-level
data are downloaded as files. The analysis reads HDF5 feature metadata via HTTP
range requests. The released objects are not the signature CSV used for fitting.
"""

import csv
import gzip
import hashlib
import json
import platform
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import fsspec
import h5py
import requests

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
TABLE = PROJECT / "panel_sources" / "GSM8558008_NEADT_12_features.tsv.gz"
ASSETS = {
    "breast": "5466ce85-7f19-49e2-9a1b-54646bb84708",
    "dlbcl": "e9b97c4b-5000-488d-b2e4-505bfb327671",
    "lung": "4ca422ae-e4ad-4884-9073-7c3799cd9143",
}
COLLECTION = "bd552f76-1f1b-43a3-b9ee-0aace57e90d6"


def string_column(group, key):
    column = group[key]
    if isinstance(column, h5py.Dataset):
        return column.asstr()[:].tolist()
    categories = column["categories"].asstr()[:]
    codes = column["codes"][:]
    if min(codes) < 0:
        raise ValueError(f"Missing categorical values in {key}")
    return categories[codes].tolist()


def audit_axis(item, lookup, duplicate_ids):
    tissue, asset = item
    url = f"https://datasets.cellxgene.cziscience.com/{asset}.h5ad"
    head = requests.head(url, timeout=25)
    head.raise_for_status()
    with fsspec.open(url, block_size=1024 * 1024) as stream, h5py.File(stream) as h5:
        var = h5["var"]
        ids = string_column(var, var.attrs["_index"])
        symbols = string_column(var, "feature_name")
        filtered = var["feature_is_filtered"][:]
        has_raw = "raw/var" in h5
        n_cells = len(h5["obs"][h5["obs"].attrs["_index"]])
    if len(ids) != len(set(ids)) or len(ids) != len(symbols):
        raise ValueError(f"Invalid feature axis for {tissue}")
    counts = Counter(lookup[gene] for gene in ids if gene in lookup)
    deposited_by_id = dict(zip(ids, symbols))
    rows = []
    for symbol, members in sorted(duplicate_ids.items()):
        for gene in members:
            rows.append(dict(tissue=tissue, stable_id=gene, symbol_2020A=symbol,
                             present=int(gene in deposited_by_id),
                             deposited_name=deposited_by_id.get(gene, ""),
                             deposited_members=counts[symbol], full_reference_members=len(members)))
    digest = hashlib.sha256()
    for gene, name, flag in zip(ids, symbols, filtered):
        digest.update(f"{gene}\t{name}\t{int(flag)}\n".encode("utf-8"))
    return dict(tissue=tissue, url=url, etag=head.headers.get("ETag"),
                axis_sha256_id_name_flag_lf=digest.hexdigest(), observations=n_cells,
                var_rows=len(ids), raw_var_present=has_raw,
                flagged_filtered=int(sum(filtered)), ids_found_2020A=sum(counts.values()),
                ids_absent_2020A=len(ids) - sum(counts.values()),
                unmatched_2020A={gene: name for gene, name in deposited_by_id.items() if gene not in lookup},
                duplicate_2020A_symbols_with_any_member=sum(counts[k] > 0 for k in duplicate_ids),
                duplicate_2020A_symbols_with_multiple_members=sum(n > 1 for n in counts.values())), rows


def main():
    with gzip.open(TABLE, "rt", encoding="utf-8") as stream:
        pairs = [(gene, symbol) for gene, symbol, kind in csv.reader(stream, delimiter="\t")
                 if kind == "Gene Expression"]
    lookup = dict(pairs)
    if len(lookup) != len(pairs):
        raise ValueError("Repeated gene ID in 2020-A comparator")
    groups = defaultdict(list)
    for gene, name in pairs:
        groups[name].append(gene)
    duplicates = {name: genes for name, genes in groups.items() if len(genes) > 1}
    with ThreadPoolExecutor(max_workers=3) as pool:
        axes = list(pool.map(lambda item: audit_axis(item, lookup, duplicates), ASSETS.items()))
    rows = [row for _, tissue_rows in axes for row in tissue_rows]
    fields = ["tissue", "stable_id", "deposited_name", "symbol_2020A", "deposited_members",
              "full_reference_members", "present"]
    with (HERE / "duplicate_group_members.tsv").open("w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(out, fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    result = dict(pmcid="PMC12069714", collection=COLLECTION, source="released Chromium var",
                  comparator=str(TABLE.relative_to(PROJECT)), comparator_rows=len(lookup),
                  comparator_duplicate_symbols=len(duplicates),
                  comparator_gzip_sha256=hashlib.sha256(TABLE.read_bytes()).hexdigest(),
                  versions=dict(python=platform.python_version(), h5py=h5py.__version__,
                                fsspec=fsspec.__version__), axes=[a for a, _ in axes],
                  model_input_identity="indeterminate")
    (HERE / "result.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
