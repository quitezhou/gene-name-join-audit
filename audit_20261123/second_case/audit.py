"""Trace published WTA/signature axes without assigning locus identity from a symbol.

Run with Python 3, h5py and numpy. Optional --refresh-atlas needs fsspec, aiohttp.
Remote archive and wrapper are pinned by SHA256. Existing four author files
must match their Git blob IDs at the pinned revision listed in result.json.
"""
import argparse
import csv
import gzip
import hashlib
import io
import json
import tarfile
from collections import defaultdict, Counter
from pathlib import Path

import h5py
import numpy as np

OUT = Path(__file__).resolve().parent
PROJ = OUT.parent.parent
SRC = OUT / "source"
ATLAS = "https://cellgeni.cog.sanger.ac.uk/covid19-nanostring/lung_atlas_integrated_RNA.cellxgene.h5ad"
ARCHIVE = "https://cellgeni.cog.sanger.ac.uk/covid19-nanostring/Bayraktar_WTA_COVID19_Atlas_2023.tar.gz"
AUTHOR_REV = "83ce1d09bca1aa003e7c951af17835d95231b9f7"
WRAPPER_REV = "dd35f3686d923294de05b55d0c59d691f56c1904"
WRAPPER_URL = f"https://raw.githubusercontent.com/BayraktarLab/cell2location/{WRAPPER_REV}/cell2location/run_c2l.py"
AUTHOR_INPUTS = {
    "UKCIC_C2L_signature_covid19.csv": ("spatial_decon/data/C2L_RegressionModel_Mean_Expression_sc_ref_covid19.csv","3572b1173d00941617f1fdd458cd850dee2307df"),
    "UKCIC_C2L_signature_healthy.csv": ("spatial_decon/data/C2L_RegressionModel_Mean_Expression_sc_ref_healthy.csv","fab746b8120fd295fd0f6b3f95cdfa24afe0a5a7"),
    "UKCIC_Run_COVID19_C2L.ipynb": ("spatial_decon/notebooks/Run_COVID19_to_COVID19_C2L.ipynb","d6bdfc289a6c72caafd8648b4fb9636b15e3d25b"),
    "UKCIC_Run_Healthy_C2L.ipynb": ("spatial_decon/notebooks/Run_Healthy_to_Healthy_C2L.ipynb","f60a86dc4fea51f83780e9e97f3574590f032ff6"),
}

def ensure_file(path, url, expected):
    if not path.exists():
        from urllib.request import urlopen
        with urlopen(url, timeout=90) as response:
            data = response.read()
        if hashlib.sha256(data).hexdigest() != expected:
            raise ValueError(f"Unexpected remote content at {url}")
        path.write_bytes(data)
    if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
        raise ValueError(f"Unexpected local content at {path}")

def sha256(data):
    return hashlib.sha256(data).hexdigest()

def dump(path, obj):
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

def audit(refresh_atlas):
    SRC.mkdir(exist_ok=True)
    archive_path = SRC / "wta.tar.gz"
    ensure_file(archive_path, ARCHIVE, "bd1972b49b2182c0d8ddae88e065a822dadd87dae11e9c5445d3f55e91e03c2b")
    wrapper_path = SRC / "run_c2l_dd35f36.py"
    ensure_file(wrapper_path, WRAPPER_URL, "58911de3851e4c5377b833a4216ddb8dbf928a6f60c304aa3c52b73fb130566a")
    wrapper = wrapper_path.read_text(encoding="utf-8")
    for expression in ("sp_ind = sp_data.var_names.isin(cell_state_df.index)",
                       "cell_state_df = cell_state_df.loc[sp_data.var_names[sp_ind], :]",
                       "X_data = X_data[:, sp_ind]"):
        assert expression in wrapper

    author_manifest = {}
    for name, (repo_path, expected) in AUTHOR_INPUTS.items():
        data = (PROJ / name).read_bytes()
        blob = hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest()
        assert blob == expected, name
        author_manifest[name] = {"url":f"https://github.com/thjimmylee/UKCIC-COVID19-paper-figures/blob/{AUTHOR_REV}/{repo_path}", "git_blob":blob, "sha256":sha256(data)}

    atlas_path = SRC / "atlas_var.json"
    if refresh_atlas or not atlas_path.exists():
        import fsspec
        with fsspec.open(ATLAS, "rb", block_size=131072) as stream, h5py.File(stream, "r") as atlas:
            group = atlas["var"]
            dump(atlas_path, {"url":ATLAS, "size_bytes":stream.size, "var_fields":list(group),
                              "has_raw": "raw" in atlas, "index_field":group.attrs["_index"],
                              "names":group[group.attrs["_index"]].asstr()[:].tolist()})
    atlas_meta = json.loads(atlas_path.read_text(encoding="utf-8"))
    atlas_names = set(atlas_meta["names"])
    ann = json.loads(gzip.decompress((PROJ / "gencode.v32.gene_metadata.json.gz").read_bytes()))
    symbol_ids = defaultdict(set)
    for gid, records in ann["gene_records"].items():
        for record in records:
            symbol_ids[record["gene_name"]].add(gid)
    result = {"pmcid":"PMC11893906", "doi":"10.1038/s41467-025-56473-x",
              "scope":"WTA host gene axes at the signature/spatial boundary. No fitted-model rerun.",
              "adjudication":"indeterminate: original index-to-locus crosswalk unavailable.",
              "wrapper_comparator_url":WRAPPER_URL,
              "wrapper_version_note":"Compatible historical source with the notebook's ignore_deprecation_exception argument. Exact installed revision unverified.",
              "archive":{"url":ARCHIVE,"sha256":sha256(archive_path.read_bytes())},
              "atlas":{"url":ATLAS,"file_bytes":atlas_meta["size_bytes"],"rows":len(atlas_meta["names"]),"unique_names":len(atlas_names),"var_fields":atlas_meta["var_fields"],"has_raw":atlas_meta["has_raw"],"index_text_sha256":sha256("\n".join(atlas_meta["names"]).encode())},
              "author_inputs":author_manifest, "annotation_comparator":{k:ann[k] for k in ["source_url","release","assembly","compressed_gtf_sha256"]},
              "arms":{}}
    trace_fields = ["arm","symbol","spatial_column_0based","signature_row_0based","atlas_var_contains_name","gencode_v32_distinct_ids","gencode_v32_candidate_ids"]
    with (OUT / "matched_feature_trace.tsv").open("w",encoding="utf-8",newline="") as output, tarfile.open(archive_path) as t:
        writer = csv.writer(output,delimiter="\t",lineterminator="\n")
        writer.writerow(trace_fields)
        for arm, stem, notebook in [("covid19","COVID-19","COVID19"),("healthy","Healthy","Healthy")]:
            member = f"AnnData_WTA_{stem}.h5ad"
            nb = json.loads((PROJ / f"UKCIC_Run_{notebook}_C2L.ipynb").read_text(encoding="utf-8"))
            code = "\n".join("".join(c["source"]) for c in nb["cells"] if c["cell_type"]=="code")
            assert member in code and "'use_raw': True" in code and "run_cell2location(" in code
            assert f"C2L_RegressionModel_Mean_Expression_sc_ref_{arm}.csv" in code
            with (PROJ / f"UKCIC_C2L_signature_{arm}.csv").open(encoding="utf-8",newline="") as sig:
                reader = csv.reader(sig)
                header = next(reader)
                ref = [r[0] for r in reader]
            sigpos = {name:i for i,name in enumerate(ref)}
            buf = t.extractfile(member).read()
            with h5py.File(io.BytesIO(buf),"r") as f:
                var = f["var"]
                names = var[var.attrs["_index"]].asstr()[:].tolist()
                rawvar = f["raw/var"]
                raw_names = rawvar[rawvar.attrs["_index"]].asstr()[:].tolist()
                assert names == raw_names and len(names)==f["raw/X"].shape[1]
                assert len(names)==len(set(names)) and len(ref)==len(sigpos)
                matched = [(i,k) for i,k in enumerate(names) if k in sigpos]
                counts = Counter()
                for i,k in matched:
                    ids = sorted(symbol_ids[k])
                    counts["absent" if not ids else "multiple" if len(ids)>1 else "single"]+=1
                    writer.writerow([arm,k,i,sigpos[k],int(k in atlas_names),len(ids),"|".join(ids)])
                details = {"spatial_file":member, "spatial_file_sha256":sha256(buf),"n_roi":f["raw/X"].shape[0],"spatial_features":len(names),
                           "spatial_var_columns":list(var.attrs["column-order"]), "raw_current_names_identical_and_ordered":True,
                           "signature_rows":len(ref), "signature_cell_states":len(header)-1,"signature_names_in_atlas":sum(k in atlas_names for k in ref),
                           "matched_symbols":len(matched),"matched_names_in_atlas":sum(k in atlas_names for _,k in matched),"matched_symbol_reference_lookup":dict(counts)}
                result["arms"][arm] = details
            with h5py.File(io.BytesIO(t.extractfile(f"AnnData_WTA_{stem}_NegProbes.h5ad").read()),"r") as f:
                neg = f["var"]; nn = neg[neg.attrs["_index"]].asstr()[:]
                details["separate_negative_control_matrix"] = {"shape":list(f["X"].shape),"unique_probe_labels":len(set(nn))}
    dump(OUT / "result.json",result)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--refresh-atlas",action="store_true")
    audit(parser.parse_args().refresh_atlas)

