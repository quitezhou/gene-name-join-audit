"""Fetch the 12 author-deposited H5s and reduce each section to exact gene totals."""
import argparse
import gzip
import hashlib
import io
import json
from pathlib import Path
import tempfile

import h5py
import numpy as np
import pandas as pd
import requests
from scipy.sparse import csc_matrix

HERE = Path(__file__).resolve().parent
BASE = "https://raw.githubusercontent.com/LieberInstitute/HumanPilot/044446d6bd8fc154aa74f7be62ec67effb1ec376/"
SAMPLES = tuple(str(i) for i in [*range(151507, 151511), *range(151669, 151677)])


def fetch(cache, url, known):
    dest = cache / url.rsplit("/", 1)[1]
    if not dest.exists():
        response = requests.get(url, timeout=90)
        response.raise_for_status()
        dest.write_bytes(response.content)
    digest = hashlib.sha256(dest.read_bytes()).hexdigest()
    if url in known and digest != known[url]:
        raise ValueError(f"File digest changed for {url}")
    return dest, dict(url=url, bytes=dest.stat().st_size, sha256=digest)


def main(cache):
    cache.mkdir(parents=True, exist_ok=True)
    manifest = HERE / "source_manifest.json"
    known = ({r["url"]: r["sha256"] for r in json.loads(manifest.read_text(encoding="utf8"))}
             if manifest.exists() else {})
    files = []
    def get(url):
        path, record = fetch(cache, url, known)
        files.append(record)
        return path
    loc = pd.read_csv(get(BASE + "AWS_File_locations.tsv"), sep="\t", dtype=str).set_index("SampleID")
    meta = pd.read_csv(get(BASE + "Analysis/visium_dlpfc_pilot_sample_metrics.tsv"), sep="\t", dtype=str).set_index("Variable")
    labels = pd.read_csv(get(BASE + "10X/barcode_level_layer_map.tsv"), sep="\t",
                         names=["barcode", "section", "layer"], dtype=str)
    assert tuple(loc.index) == SAMPLES and not labels.duplicated(["barcode", "section"]).any()
    base_axis, stats, totals, layer_data = None, [], {}, {}
    for section in SAMPLES:
        path = get(loc.at[section, "h5_filtered"])
        with h5py.File(path) as h5:
            m = h5["matrix"]
            types = m["features/feature_type"].asstr()[:]
            selected = np.where(types == "Gene Expression")[0]
            ids = list(m["features/id"].asstr()[:][selected])
            names = list(m["features/name"].asstr()[:][selected])
            assert len(set(ids)) == len(ids)
            axis = list(zip(ids,names))
            if base_axis is None:
                base_axis = axis
            elif axis != base_axis:
                raise ValueError(f"Different complete axis in {section}")
            mat = csc_matrix((m["data"][:].astype(np.int64), m["indices"][:], m["indptr"][:]),
                             shape=tuple(m["shape"][:]))[selected].tocsc()
            barcodes = m["barcodes"].asstr()[:]
            assert len(set(barcodes)) == mat.shape[1] == int(meta.at["Estimated.Number.of.Spots",section].replace(",",""))
        reference = labels[labels.section == section].set_index("barcode")["layer"]
        assert set(reference.index) <= set(barcodes)
        layer = reference.reindex(barcodes)
        totals[section] = np.asarray(mat.sum(axis=1)).ravel()
        counts = {}
        for key in ["L3","WM"]:
            mask = layer.eq(key).to_numpy()
            counts["n_" + key] = int(mask.sum())
            assert mask.any()
            layer_data[f"{section}_{key}"] = np.asarray((mat[:,mask] > 0).sum(axis=1)).ravel()
        stats.append(dict(section=section, donor=meta.at["Brain.Number",section], spots=len(barcodes),
                          total_umi=int(totals[section].sum()), annotated_spots=len(reference), **counts))
    # Gene sums are needed for all rows to reconstruct keys after filtering.
    # Layer counts are needed only for rows whose source symbol occurs more than once.
    frame = pd.DataFrame(base_axis, columns=["id","symbol"])
    compact = pd.concat([frame,pd.DataFrame(totals)],axis=1)
    layer_frame = pd.concat([frame,pd.DataFrame(layer_data)],axis=1)
    for name, out in [("gene_totals",compact),
                      ("duplicate_family_detections",layer_frame[frame.symbol.duplicated(keep=False)])]:
        with (HERE / f"{name}.tsv.gz").open("wb") as target:
            with gzip.GzipFile(filename="", fileobj=target, mode="wb", mtime=0) as compressed:
                with io.TextIOWrapper(compressed, encoding="utf8", newline="") as handle:
                    out.to_csv(handle, sep="\t",index=False,lineterminator="\n")
    pd.DataFrame(stats).to_csv(HERE/"sections.tsv",sep="\t",index=False,lineterminator="\n")
    manifest.write_text(json.dumps(files,indent=2)+"\n",encoding="utf8")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--cache",type=Path,default=Path(tempfile.gettempdir())/"gb_ext_dlpfc")
    main(parser.parse_args().cache)
