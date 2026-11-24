"""Construct cellwise witnesses and enumerate the finite source-stratum extreme points."""
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent


def main():
    x = pd.read_csv(HERE / "cell_detection_trace.tsv.gz", sep="\t")
    summary = json.loads((HERE / "contrast_bounds.json").read_text())
    df = pd.read_csv(HERE.parent/"hwat_model/source/ly_rows.tsv.gz",sep="\t")
    assert df.cell.equals(x.cell)
    assert np.array_equal(df["counts_POLR2J3"] > 0, x.bare)
    assert np.array_equal(df["counts_POLR2J3.1"] > 0, x.printed)
    assert np.array_equal(df.cluster == 3, x.in_cluster)
    for arm, fixed in [("published_followup_baseline", x.before == 1), ("all_verified", x.resolved >= 0)]:
        y = summary[arm]
        inside = x.in_cluster.astype(bool)
        def diff(v):
            return 100 * (v[inside].sum()/inside.sum()-v[~inside].sum()/(~inside).sum())
        def candidate(w):
            return np.where(fixed, x.resolved, w)
        local_lo = np.minimum(x.bare, x.printed)
        local_hi = np.maximum(x.bare, x.printed)
        for low, high, k in [(0,1,"unrestricted_binary_interval_pp"),
                              (local_lo,local_hi,"two_row_cell_interval_pp")]:
            bounds = [diff(candidate(np.where(inside,low,high))),diff(candidate(np.where(inside,high,low)))]
            np.testing.assert_allclose(bounds, y[k], atol=1e-9, rtol=0)
        ni,no = int(inside.sum()),int((~inside).sum())
        fixed_total = (x.loc[inside&fixed,"resolved"].sum()*no
                       -x.loc[~inside&fixed,"resolved"].sum()*ni)
        alternatives=[]
        for _, block in x.groupby("source_stratum"):
            known=block.loc[fixed[block.index]]
            unknown=block.loc[~fixed[block.index]]
            valid=[key for key in ["printed","bare"] if known[key].eq(known.resolved).all()]
            possible=set(int(unknown.loc[unknown.in_cluster==1,key].sum()*no
                           -unknown.loc[unknown.in_cluster==0,key].sum()*ni) for key in valid)
            assert possible
            alternatives.append(possible)
        values=[(fixed_total+sum(choices))*100/(ni*no) for choices in itertools.product(*alternatives)]
        np.testing.assert_allclose([min(values),max(values)],y["two_row_source_stratum_interval_pp"],atol=1e-9,rtol=0)
        assert (x.loc[fixed,"printed"] != x.loc[fixed,"resolved"]).sum() == y["observed_discordant_cells"]["all"]
    assert x.before.sum() == 1375 and (x.resolved>=0).sum() == 1476
    assert int(((x.before==0) & (x.resolved>=0) & (x.printed!=x.resolved)).sum())==7
    print("Reconciled 10,394 source-object rows, both arm counts, six bounds and 273 discordant detections.")


if __name__ == "__main__":
    main()
