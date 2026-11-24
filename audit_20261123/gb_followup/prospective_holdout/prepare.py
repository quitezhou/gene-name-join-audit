"""Build public per-feature aggregates from the original GEO matrices and singlet calls."""
import argparse
import csv
import gzip
import io
import json
from hashlib import sha256
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.io import mmread

HERE = Path(__file__).resolve().parent
CONDITIONS = {"ctrl": ("GSM2560248", "2.1"), "stim": ("GSM2560249", "2.2")}


def main(cache):
    files = json.loads((HERE / "sources.json").read_text())
    for source in files:
        path = cache / source["name"]
        if sha256(path.read_bytes()).hexdigest() != source["sha256"]:
            raise ValueError(f"Incorrect source bytes: {path}")
    axis = pd.read_csv(cache / "GSE96583_batch2.genes.tsv.gz", header=None,
                       names=["id", "symbol"], sep="\t", keep_default_na=False)
    assert axis.id.is_unique and not (axis == "").any().any()
    metadata = pd.read_csv(cache / "GSE96583_batch2.total.tsne.df.tsv.gz",
                           index_col=0, sep="\t", dtype={"ind": str})
    samples = []
    for condition, (gsm, mtx) in CONDITIONS.items():
        barcodes = pd.read_csv(cache / f"{gsm}_barcodes.tsv.gz", header=None)[0]
        subset = metadata.loc[metadata.stim == condition].copy()
        # In the pooled metadata 313 stim barcodes have "-11" instead of "-1".
        # Within a condition their DNA prefix is unique; explicitly assert the
        # DNA sequence for every matrix column in its supplied order.
        source_dna = barcodes.str.split("-").str[0].to_list()
        meta_dna = subset.index.str.split("-").str[0].to_list()
        assert len(set(source_dna)) == len(source_dna) and source_dna == meta_dna
        with gzip.open(cache / f"{gsm}_{mtx}.mtx.gz", "rb") as fin:
            mat = mmread(fin).tocsc()
        assert mat.shape == (len(axis), len(subset))
        assert np.all(mat.data > 0) and np.all(mat.data == np.floor(mat.data))
        mat = mat.astype(np.int64)
        mat.sum_duplicates()
        for donor in sorted(subset.ind.dropna().unique(), key=int):
            selected = ((subset.ind == donor) & (subset.multiplets == "singlet")).to_numpy()
            counts = mat[:, selected].tocsr()
            name = f"{condition}_{donor}"
            axis[name + "_cells"] = np.diff(counts.indptr)
            axis[name + "_umis"] = np.asarray(counts.sum(axis=1)).ravel()
            samples.append((name, condition, donor, int(selected.sum()), int(counts.sum())))
    with (HERE / "samples.tsv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(("sample", "condition", "donor", "singlets", "UMIs"))
        writer.writerows(samples)
    with (HERE / "axis_counts.tsv.gz").open("wb") as raw:
        with gzip.GzipFile(fileobj=raw, filename="", mode="wb", mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as stream:
                axis.to_csv(stream, sep="\t", lineterminator="\n", index=False)


if __name__ == "__main__":
    arg = argparse.ArgumentParser(description=__doc__)
    arg.add_argument("cache", type=Path)
    main(arg.parse_args().cache)
