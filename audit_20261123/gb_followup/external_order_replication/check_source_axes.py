"""Cross-check raw H5 totals and keys using AnnData; export source-layer UMI sums."""
import argparse
import json
from pathlib import Path
import tempfile
import sys

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from extract import BASE, HERE, fetch
sys.path.insert(0, str(HERE.parent / "systematic_mechanism"))
from audit_join import assign_keys


def main(cache):
    from scipy.sparse import csc_matrix
    manifest = json.loads((HERE/"source_manifest.json").read_text())
    files = {r["url"]: r["sha256"] for r in manifest}
    old = HERE/"model_sources.json"
    if old.exists():
        files.update({r["url"]:r["sha256"] for r in json.loads(old.read_text())})
    results = json.loads((HERE/"summary.json").read_text())
    names = results["strict_keys"]
    model_sources = []
    for file in ["rda/sce_layer.Rdata", "rda/eb_contrasts.Rdata", "layer_specificity.R"]:
        _,source = fetch(cache, BASE+"Analysis/Layer_Guesses/"+file,files)
        model_sources.append(source)
    old.write_text(json.dumps(model_sources,indent=2)+"\n",encoding="utf8")
    labels_file,_ = fetch(cache,BASE+"10X/barcode_level_layer_map.tsv",files)
    labels = pd.read_csv(labels_file,sep="\t",names=["barcode","section","layer"],dtype=str)
    ref = pd.read_csv(HERE/"gene_totals.tsv.gz",sep="\t",keep_default_na=False)
    sig = pd.read_csv(HERE/"layer_contrasts.tsv",sep="\t",keep_default_na=False,dtype={"section":str})
    output = []
    for r in manifest:
        if not r["url"].endswith(".h5"):
            continue
        path,_ = fetch(cache,r["url"],files)
        section = path.name.split("_")[0]
        with h5py.File(path) as h5:
            m = h5["matrix"]
            ids = m["features/id"].asstr()[:]
            symbols = m["features/name"].asstr()[:]
            assert set(m["features/feature_type"].asstr()[:]) == {"Gene Expression"}
            bar = m["barcodes"].asstr()[:]
            counts = m["data"][:]
            idx = m["indices"][:]
            all_sums = np.bincount(idx,weights=counts,minlength=len(ids)).astype(np.int64)
            assert (all_sums == ref[section]).all() and (ids == ref.id).all() and (symbols == ref.symbol).all()
            sparse = csc_matrix((counts,idx,m["indptr"][:]),shape=tuple(m["shape"][:]))
        masks = {"protected":np.ones(len(ids),dtype=bool),"filtered":all_sums>0}
        obs_layer = labels[labels.section==section].set_index("barcode").layer.reindex(bar)
        for mode,mask in masks.items():
            obj = ad.AnnData(shape=(0,int(mask.sum())))
            obj.var_names = symbols[mask]
            obj.var_names_make_unique()
            assert list(obj.var_names) == list(assign_keys(symbols[mask], "-"))
            index = dict(zip(obj.var_names,np.flatnonzero(mask)))
            # All members' keys under full-axis naming, including suffixed records.
            if mode == "protected":
                anchors = {key:ids[i] for key,i in index.items() if symbols[i] in names}
            for key,anchor in anchors.items():
                original_row = index.get(key)
                source_id = ids[original_row] if original_row is not None else ""
                srows = sig[(sig.section==section)&(sig.arm==mode)&(sig.key==key)]
                vector = (sparse[original_row].toarray()[0] if original_row is not None
                          else np.zeros(sparse.shape[1],dtype=np.int64))
                for level in sorted(obs_layer.dropna().unique()):
                    layer_mask = obs_layer.eq(level).to_numpy()
                    umi = int(vector[layer_mask].sum())
                    n = int((vector[layer_mask]>0).sum())
                    if len(srows) and level in ("L3","WM"):
                        assert source_id == srows.iloc[0].id and n == srows.iloc[0]["detected_"+level]
                    layer_name = level.replace("L","Layer") if level.startswith("L") else level
                    output.append(dict(section_layer=section+"_"+layer_name,key=key,anchor_id=anchor,
                                       source_id=source_id,arm=mode,umi=umi))
    pd.DataFrame(output).to_csv(HERE/"family_layer_counts.tsv",sep="\t",index=False)
    print("Independent H5 gene totals and AnnData marker projections agree across all 12 sections.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--cache",type=Path,default=Path(tempfile.gettempdir())/"gb_ext_dlpfc")
    main(parser.parse_args().cache)
