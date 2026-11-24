"""Check every reconstructed score via the CellTypist 1.2.0 public annotate API.

With celltypist==1.2.0 installed call python check_celltypist.py.
Alternatively pass --wheel PATH to its official 1.2.0 distribution wheel.
Run projection.py first to prepare expression arrays and independent fit scores.
"""
import argparse
import logging
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import mmread

from projection import COMMON, HERE, MODEL_BLOBS, OBJECTS, read_model


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--wheel", type=Path)
    args = parser.parse_args()
    if args.wheel:
        sys.path.insert(0, str(args.wheel.resolve()))
    import anndata
    import celltypist
    from celltypist.models import Model
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler

    assert celltypist.__version__ == "1.2.0"
    logging.getLogger("celltypist").setLevel(logging.ERROR)
    cache = Path(tempfile.gettempdir()) / "massier_endpoint_projection"
    models = {}
    for name in MODEL_BLOBS:
        lr, sc, _ = read_model(name)
        regression, scaler = LogisticRegression(), StandardScaler()
        regression.__dict__.update(vars(lr))
        scaler.__dict__.update(vars(sc))
        models[name] = Model(regression, scaler, {})
    rows = []
    for obj in OBJECTS:
        keys = (cache/f"{obj}_keys.txt").read_text().splitlines()
        counts = mmread(cache/f"{obj}.mtx").tocsr()
        meta = pd.read_csv(cache/f"{obj}_obs.tsv", sep="\t")
        for source, obs in meta.groupby("source"):
            for name, model in models.items():
                file = cache/f"{obj}_{source}_{name}_scores.npz"
                expected = np.load(file)
                selected = sorted(set(keys) & set(model.classifier.features))
                # Shuffling gene order relative to the projection tests name-based alignment.
                outcol = [keys.index(s) for s in selected]
                labels = selected + ["__other_umi_for_size_check__"]
                for arm in ["baseline", "aligned"]:
                    observed = counts[outcol][:, obs.index].T.toarray()
                    if arm == "aligned":
                        for key in set(selected) & COMMON.keys():
                            target = COMMON[key][2][source]
                            observed[:,selected.index(key)] = (counts[keys.index(target),obs.index].toarray()[0]
                                                              if target in keys else 0.)
                    extra = obs.total.to_numpy() - observed.sum(axis=1)
                    assert min(extra) >= -1e-7
                    raw = np.column_stack([observed, np.maximum(extra, 0.)])
                    x = np.log1p(10000.0 * raw / obs.total.to_numpy()[:,None])
                    adata = anndata.AnnData(X=x, obs=pd.DataFrame(index=obs.barcode),
                                            var=pd.DataFrame(index=labels))
                    observed_result = celltypist.annotate(adata, model=model, majority_voting=False)
                    df = observed_result.decision_matrix.loc[obs.barcode, model.classifier.classes_]
                    err = np.abs(df.to_numpy() - expected[arm])
                    labs = observed_result.predicted_labels.loc[obs.barcode, "predicted_labels"].astype(str)
                    truth = model.classifier.classes_[expected[arm].argmax(axis=1)]
                    assert np.allclose(df.to_numpy(), expected[arm], atol=1e-9, rtol=1e-10)
                    assert np.all(labs==truth)
                    rows.append(dict(object=obj, source=source, model=name, arm=arm,
                                     n_cells=len(obs), n_classes=len(model.classifier.classes_),
                                     max_absolute_score_difference=err.max(), label_disagreements=int((labs!=truth).sum()),
                                     celltypist=celltypist.__version__))
    pd.DataFrame(rows).to_csv(HERE/"official_api_verification.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
