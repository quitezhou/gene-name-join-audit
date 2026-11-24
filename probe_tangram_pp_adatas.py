"""Execute unmodified pp_adatas AST from a pinned official Tangram source file.

This is a source-only function probe with real AnnData/NumPy/SciPy objects.
Scanpy package import stalled, so its unmodified filter_genes body and the
required in-memory helpers are extracted from the installed source by AST.
This is not an installed Tangram/Scanpy package integration test.
"""
import ast
import hashlib
import importlib.metadata
import json
import logging
from functools import singledispatch
from pathlib import Path
from types import SimpleNamespace
from urllib.request import urlopen

import anndata
import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

COMMIT = "4c68995a418f41dc8caef567598c4d9b47781a13"
URL = f"https://raw.githubusercontent.com/broadinstitute/Tangram/{COMMIT}/tangram/mapping_utils.py"


def source_only_scanpy():
    package = Path(importlib.metadata.distribution("scanpy").locate_file("scanpy"))
    files = {"filter_genes": package / "preprocessing/_simple.py",
             "materialize_as_ndarray": package / "preprocessing/_distributed.py",
             "axis_sum": package / "_utils/__init__.py",
             "is_backed_type": package / "_utils/__init__.py",
             "raise_not_implemented_error_if_backed_type": package / "_utils/__init__.py"}
    scope = {"np": np, "AnnData": anndata.AnnData, "h5py": h5py,
             "SparseDataset": anndata.abc.CSRDataset | anndata.abc.CSCDataset,
             "DaskArray": type("UnusedDaskArrayPlaceholder", (), {}),
             "singledispatch": singledispatch,
             "old_positionals": lambda *args: lambda fn: fn,
             "logg": logging}
    hashes = {}
    for name, path in files.items():
        raw = path.read_bytes()
        hashes[str(path)] = hashlib.sha256(raw).hexdigest()
        nodes = [n for n in ast.parse(raw.decode("utf-8")).body
                 if isinstance(n, ast.FunctionDef) and n.name == name]
        func = nodes[-1]  # Last is implementation, earlier nodes are overloads.
        future = ast.ImportFrom(module="__future__", names=[ast.alias(name="annotations")], level=0)
        module = ast.fix_missing_locations(ast.Module(body=[future, func], type_ignores=[]))
        exec(compile(module, str(path), "exec"), scope)
    return SimpleNamespace(pp=SimpleNamespace(filter_genes=scope["filter_genes"])), hashes


def make_data(stable_names=False, spatial_zero=False):
    ids = ["ENSG_TEST_1", "ENSG_TEST_2", "ENSG_TEST_3"]
    symbols = ["DUP", "DUP", "CONTROL"]
    var = pd.DataFrame({"stable_id": ids, "original_symbol": symbols},
                       index=ids if stable_names else symbols)
    ref = anndata.AnnData(csr_matrix([[2, 3, 5], [4, 5, 6]], dtype=float), var=var.copy())
    spatial = anndata.AnnData(csr_matrix([[0 if spatial_zero else 1, 7, 5],
                                        [0 if spatial_zero else 2, 8, 6]], dtype=float),
                             var=var.copy())
    return ref, spatial


def run():
    source_path = Path(__file__).with_name("tangram_mapping_utils_pinned.py")
    raw = source_path.read_bytes() if source_path.exists() else urlopen(URL, timeout=60).read()
    tree = ast.parse(raw.decode("utf-8"))
    func = next(node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == "pp_adatas")
    sc, scanpy_hashes = source_only_scanpy()
    scope = {"np": np, "pd": pd, "sc": sc, "logging": logging}
    exec(compile(ast.Module(body=[func], type_ignores=[]), URL, "exec"), scope)
    rows = []
    for name, stable, zero, lowercase in [
        ("same_feature_table_all_expressed", False, False, True),
        ("same_feature_table_differential_zero_filter", False, True, True),
        ("differential_zero_filter_without_lowercase", False, True, False),
        ("stable_id_index_control", True, True, True),
    ]:
        ref, spatial = make_data(stable_names=stable, spatial_zero=zero)
        before_identical = ref.var.equals(spatial.var)
        scope["pp_adatas"](ref, spatial, gene_to_lowercase=lowercase)
        selected = list(ref.uns["training_genes"])
        pairs = [{"training_key": key,
                  "reference_id": str(ref.var.loc[key, "stable_id"]),
                  "spatial_id": str(spatial.var.loc[key, "stable_id"])}
                 for key in selected]
        mismatches = [p for p in pairs if p["reference_id"] != p["spatial_id"]]
        rows.append({"case": name, "input_feature_tables_identical": before_identical,
                     "reference_keys_after": list(ref.var_names),
                     "spatial_keys_after": list(spatial.var_names),
                     "training_keys_agree": selected == list(spatial.uns["training_genes"]),
                     "pairs": pairs, "mismatch_count": len(mismatches), "mismatches": mismatches})
    assert [x["mismatch_count"] for x in rows] == [0, 1, 1, 0]
    result = {"scope": "source-only AST execution of unmodified Tangram pp_adatas and installed Scanpy filter_genes/helper bodies; actual AnnData/NumPy/SciPy; no Tangram optimizer run; not an installed-package integration test",
              "import_limitation": "Full Scanpy import stalled; a bounded faulthandler attempt gave an incomplete traceback. Original processes were stopped.",
              "scanpy_source_sha256": scanpy_hashes,
              "adapter_scope": "Keyword-only calls; old_positionals compatibility decorator replaced by identity, DaskArray by unused placeholder, Python logging used. All test matrices are in-memory SciPy CSR; no Dask/backed arrays tested.",
              "commit": COMMIT, "source_url": URL,
              "source_sha256": hashlib.sha256(raw).hexdigest(),
              "versions": {p: importlib.metadata.version(p) for p in ["scanpy", "anndata", "numpy", "pandas", "scipy"]},
              "cases": rows}
    output = Path(__file__).with_name("tangram_pp_adatas_probe_result.json")
    output.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    run()
