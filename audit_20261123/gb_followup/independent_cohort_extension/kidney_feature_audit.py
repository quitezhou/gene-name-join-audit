"""Audit first-record retention on the published CKD inputs. No cell2location fit is inferred."""
import argparse
from collections import Counter
from hashlib import md5, sha256
import io
import json
from pathlib import Path
import tarfile
import zipfile

import h5py
import numpy as np
import pandas as pd
import requests
from scipy.io import mmread
import gzip


SOURCES = {
    "spatial.zip": ("https://ndownloader.figshare.com/files/42173859", "e8a3768d6bd05db266a343c3c240357d"),
    "code.zip": ("https://ndownloader.figshare.com/files/42171600", "e46b3aa81f3c06e81ade76c64a3b6b8a"),
    "Human_CD10negative.tar.gz": ("https://zenodo.org/api/records/4059315/files/Human_CD10negative.tar.gz/content",
                                 "4af0d2909c2fb1053f3e772eaebedefe"),
}
HERE = Path(__file__).resolve().parent


def obtain(cache, name):
    url, expected = SOURCES[name]
    path = cache / name
    if not path.exists():
        with requests.get(url, stream=True, timeout=120) as resp, path.with_suffix(".part").open("wb") as out:
            resp.raise_for_status()
            for part in resp.iter_content(1024*1024):
                out.write(part)
        path.with_suffix(".part").replace(path)
    check = md5()
    with path.open("rb") as src:
        for part in iter(lambda: src.read(4*1024*1024), b""):
            check.update(part)
    assert check.hexdigest() == expected
    return path


def read_reference(path):
    with tarfile.open(path, "r:gz") as src:
        for mem in src:
            if mem.name.endswith("colData.txt"):
                cells = pd.read_csv(src.extractfile(mem))
            elif mem.name.endswith("rowData.txt"):
                feat = pd.read_csv(src.extractfile(mem), keep_default_na=False)
            elif mem.name.endswith(".mtx"):
                f = src.extractfile(mem)
                assert f.readline().strip() == b"%%MatrixMarket matrix coordinate integer general"
                nr, nc, nonzero = map(int, f.readline().split())
                assert (nr, nc) == (len(feat), len(cells))
                totals = np.zeros(nr, dtype=np.int64)
                detected = np.zeros(nr, dtype=np.int64)
                prev_col = prev_row = -1
                n = 0
                for chunk in pd.read_csv(f, sep=" ", header=None, names=["row","col","val"], chunksize=1500000):
                    row, col, val = (chunk[k].to_numpy() for k in ("row", "col", "val"))
                    # The deposited matrix is column-major and unique within each column.
                    assert (col[0], row[0]) > (prev_col, prev_row)
                    assert np.all((col[1:] > col[:-1]) | ((col[1:] == col[:-1]) & (row[1:] > row[:-1])))
                    assert (val > 0).all()
                    prev_col, prev_row = col[-1], row[-1]
                    detected += np.bincount(row - 1, minlength=nr)
                    totals += np.bincount(row - 1, weights=val, minlength=nr).astype(np.int64)
                    n += len(row)
                assert nonzero == n
                break
    feat.columns = ["symbol", "stable_id"]
    assert feat.stable_id.is_unique
    feat["cells"] = detected
    feat["umi"] = totals
    return feat, nc


def read_spatial(path):
    sections = []
    fingerprints = {}
    with zipfile.ZipFile(path) as z:
        for h5 in sorted(n for n in z.namelist() if n.endswith(".h5")):
            prefix = h5.rsplit("/",1)[0] + "/"
            sample = h5.split("/")[0]
            feat_bytes = gzip.decompress(z.read(prefix + "features.tsv.gz"))
            fingerprints[sample] = sha256(feat_bytes).hexdigest()
            tab = pd.read_csv(io.BytesIO(feat_bytes), sep="\t", header=None, names=["stable_id", "symbol", "type"])
            with h5py.File(io.BytesIO(z.read(h5))) as f:
                h = f["matrix"]
                assert tab.stable_id.tolist() == h["features/id"].asstr()[:].tolist()
                assert tab.symbol.tolist() == h["features/name"].asstr()[:].tolist()
                calc = np.bincount(h["indices"][:], weights=h["data"][:], minlength=len(tab)).astype(np.int64)
                assert tab.stable_id.is_unique
            with gzip.open(io.BytesIO(z.read(prefix + "matrix.mtx.gz"))) as stream:
                mex = mmread(stream).tocsr()
            assert np.array_equal(np.asarray(mex.sum(axis=1)).ravel(), calc)
            tab["umi"] = calc
            tab["spots"] = np.diff(mex.indptr)
            sections.append(tab)
    assert len(sections) == 2 and len(set(fingerprints.values())) == 1
    return sections, fingerprints


def compare(cache):
    paths = {k: obtain(cache, k) for k in SOURCES}
    with zipfile.ZipFile(paths["code.zip"]) as z:
        code = z.read("codes/Fig1/Fig1_process_data.R").decode()
        for expr in ["duplicated(gene.annot[,2])", "mat <- mat[-to.rm,]", "rownames(mat) <- gene.annot[,2]",
                     "CreateSeuratObject(mat,min.cells = 50)"]:
            assert expr in code
    ref, n_cells = read_reference(paths["Human_CD10negative.tar.gz"])
    spatial, fingerprints = read_spatial(paths["spatial.zip"])
    ref.to_csv(HERE / "kidney_reference_axis.tsv.gz", sep="\t", index=False, compression={"method":"gzip","mtime":0})
    spatial_export = spatial[0][["symbol","stable_id"]].copy()
    for i, s in enumerate(spatial, 1):
        spatial_export[f"umi{i}"] = s.umi
        spatial_export[f"spots{i}"] = s.spots
    spatial_export.to_csv(HERE / "kidney_spatial_axis.tsv.gz", sep="\t", index=False, compression={"method":"gzip","mtime":0})
    sc_by_id = ref.set_index("stable_id")
    by_symbol = ref.drop_duplicates("symbol").set_index("symbol")
    multiplicity = ref.symbol.value_counts()
    # Use the bare-name row; the reference script retains only the first row per symbol.
    v = [s.drop_duplicates("symbol").set_index("symbol") for s in spatial]
    assert v[0].index.is_unique and v[1].index.is_unique
    common = sorted(set(by_symbol.index) & set(v[0].index))
    lost = n_cells - 49226
    output = []
    for gene in common:
        a, b = by_symbol.loc[gene], v[0].loc[gene]
        if a.stable_id == b.stable_id:
            continue
        same = sc_by_id.loc[b.stable_id] if b.stable_id in sc_by_id.index else None
        rec = {"symbol": gene, "ref_symbol_rows": int(multiplicity[gene]),
               "spatial_symbol_rows": int((spatial[0].symbol == gene).sum()), "ref_first_id": a.stable_id,
               "spatial_id": b.stable_id, "reference_first_umi": int(a.umi), "reference_first_detected": int(a.cells),
               "reference_same_spatial_id_symbol": same.symbol if same is not None else "",
               "reference_same_spatial_id_umi": int(same.umi) if same is not None else None,
               "reference_same_spatial_id_detected": int(same.cells) if same is not None else None}
        for i, t in enumerate(v,1):
            rec[f"spatial{i}_umi"] = int(t.at[gene, "umi"])
            rec[f"spatial{i}_spots"] = int(t.at[gene, "spots"])
        output.append(rec)
    out = pd.DataFrame(output)
    out.to_csv(HERE / "kidney_shared_name_mismatches.tsv", sep="\t", index=False, float_format="%.0f")
    dropped = ref.loc[ref.symbol.duplicated()].copy()
    dropped.to_csv(HERE / "kidney_second_and_later_records.tsv", sep="\t", index=False)
    freq = Counter(ref.symbol)
    result = dict(reference_records=len(ref), reference_cells=n_cells, published_n=49226, difference_cells=lost,
                  spatial_records=len(spatial[0]), spatial_unique_symbols=len(v[0]),
                  spatial_axis_sha256=fingerprints, shared_names=len(common), first_name_different_id=len(out),
                  duplicate_related=int(((out.ref_symbol_rows>1)|(out.spatial_symbol_rows>1)).sum()),
                  duplicated_reference_families=sum(x>1 for x in freq.values()),
                  omitted_reference_records=len(dropped), omitted_umi=int(dropped.umi.sum()),
                  detected_in_50_cells_among_omitted=int((dropped.cells>=50).sum()),
                  detected_in_50_even_if_all_lost_cells_expressed=int((dropped.cells>=50+lost).sum()),
                  source_matrices_md5={k:v[1] for k,v in SOURCES.items()})
    (HERE/"kidney_summary.json").write_text(json.dumps(result, indent=2)+"\n", encoding="utf8")
    print(json.dumps(result))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache", type=Path, required=True)
    cache = parser.parse_args().cache
    cache.mkdir(exist_ok=True, parents=True)
    compare(cache)
