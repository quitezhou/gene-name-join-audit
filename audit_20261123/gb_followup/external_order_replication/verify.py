"""Offline verification of the DLPFC report and all 22,331 BH decisions."""
from pathlib import Path
import tempfile

import numpy as np
import pandas as pd

from replay import HERE, run


def bh(values):
    p = np.asarray(values)
    rank = np.argsort(p)
    adjusted = np.empty_like(p)
    adjusted[rank] = np.minimum(1, np.minimum.accumulate((len(p)*p[rank]/np.arange(1,len(p)+1))[::-1])[::-1])
    return adjusted


def main():
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp)
        run(out)
        for name in ["summary.json","pairs.tsv","events.tsv","layer_contrasts.tsv","layer_donors.tsv"]:
            # Normalize platform line endings for the numeric TSV/JSON comparison.
            assert (out/name).read_text(encoding="utf8") == (HERE/name).read_text(encoding="utf8"),name
    fit = pd.read_csv(HERE/"full_layer_fit.tsv.gz",sep="\t",index_col="id",float_precision="round_trip")
    target = pd.read_csv(HERE/"layer_model.tsv",sep="\t").set_index(["id","arm"])
    check = pd.read_csv(HERE/"layer_model_check.tsv",sep="\t").iloc[0]
    assert len(fit) == check.genes == 22331 and fit.index.is_unique
    q = {arm:bh(fit[arm+"_p"]) for arm in ["paper","protected","counterfactual"]}
    assert all(np.isfinite(scores).all() for scores in q.values())
    assert (np.abs(fit.paper_beta-fit.protected_beta).max() < 2e-13
            and np.abs(fit.paper_p-fit.protected_p).max() < 2e-13)
    for arm in q:
        assert (q[arm]<.05).sum() == 8500
        assert np.array_equal(q[arm]<.05,q["paper"]<.05)
    for (gene,arm),row in target.iterrows():
        key = {"published":"paper","protected_replay":"protected","filter_then_name":"counterfactual"}[arm]
        assert np.isclose(row.p_BH,q[key][fit.index.get_loc(gene)],rtol=1e-10)
        assert np.isclose(row.coefficient,fit.at[gene,key+"_beta"],rtol=1e-12)
        assert np.isclose(row.p,fit.at[gene,key+"_p"],rtol=1e-12)
    assert check.checked_source_counts == 152 and check.max_error_counts == check.changed_calls == 0
    inputs = pd.read_csv(HERE/"archive_count_check.tsv",sep="\t")
    raw = pd.read_csv(HERE/"family_layer_counts.tsv",sep="\t")
    protected = raw[raw.arm=="protected"].rename(columns={"anchor_id":"id"})
    rows = inputs.merge(protected,on=["id","section_layer"],validate="one_to_one")
    assert len(rows) == len(inputs) == len(protected) == 152
    assert (rows.H5_sum == rows.archive_count).all() and (rows.umi == rows.H5_sum).all()
    print("All 22,331 WM–Layer 3 decisions reconciled; 8,500 pass the threshold in all three arms.")


if __name__=="__main__":
    main()
