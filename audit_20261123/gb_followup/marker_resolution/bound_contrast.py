"""Evaluate attainable full-object detection differences given unresolved target identities."""
import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent


def contrasts(df, known, stratum):
    cluster = df.in_cluster == 1
    per_detection = pd.Series(np.where(cluster, 100 / cluster.sum(), -100 / (~cluster).sum()), index=df.index)
    fixed = float((df.loc[known, "resolved"] * per_detection[known]).sum())
    unknown = df.loc[~known]
    weights = per_detection[~known]
    def interval(lo, hi):
        return [fixed + float(np.minimum(lo*weights, hi*weights).sum()),
                fixed + float(np.maximum(lo*weights, hi*weights).sum())]
    binary = interval(0, 1)
    two_row = interval(unknown.bare, unknown.printed)
    group_report = []
    for group, rows in df.groupby(stratum, dropna=False, sort=True):
        group_verified = rows.loc[known[rows.index] & rows.bare.ne(rows.printed)]
        possible = [key for key in ["bare", "printed"] if (group_verified[key] == group_verified.resolved).all()]
        # A mixed-source label may exhibit two verified bindings; treating that
        # label as homogeneous then conflicts with observed source records.
        if not possible:
            raise ValueError(f"One binding contradicts verified cells in {group}")
        un = rows.loc[~known[rows.index]]
        values = [(un[key]*per_detection[un.index]).sum() for key in possible]
        group_report.append(dict(stratum=group, unresolved_n=len(un), verified_n=int(known[rows.index].sum()),
                                 possible="|".join(possible), bare_pp=float((un.bare*per_detection[un.index]).sum()),
                                 printed_pp=float((un.printed*per_detection[un.index]).sum()),
                                 lower=min(values), upper=max(values)))
    g = pd.DataFrame(group_report)
    grouped = [fixed + g.lower.sum(), fixed + g.upper.sum()]
    for bounds in [binary, two_row, grouped]:
        assert bounds[0] <= bounds[1] + 1e-9
    assert binary[0]-1e-9 <= two_row[0] <= grouped[0]+1e-9
    assert grouped[1]-1e-9 <= two_row[1] <= binary[1]+1e-9
    repaired = df.printed.copy()
    repaired.loc[known] = df.loc[known, "resolved"]
    fixed_changes = known & df.printed.ne(df.resolved)
    could_change = ~known & df.bare.ne(df.printed)
    def result(v):
        ni = int(v[cluster].sum()); no = int(v[~cluster].sum())
        return {"in_positive":ni, "in_total":int(cluster.sum()), "out_positive":no,
                "out_total":int((~cluster).sum()), "contrast_pp":float((v*per_detection).sum())}
    return dict(n_verified=int(known.sum()), cluster_verified=int((known & cluster).sum()),
                released=result(df.printed), partial=result(repaired),
                observed_discordant_cells={"all":int(fixed_changes.sum()),
                                           "in_cluster":int((fixed_changes & cluster).sum()),
                                           "elsewhere":int((fixed_changes & ~cluster).sum())},
                possible_discordant_cells_under_two_row_model=[
                    int(fixed_changes.sum()), int((fixed_changes | could_change).sum())],
                unrestricted_binary_interval_pp=binary, two_row_cell_interval_pp=two_row,
                two_row_source_stratum_interval_pp=grouped), g


def main():
    df = pd.read_csv(HERE/"cell_detection_trace.tsv.gz", sep="\t")
    assert len(df) == len(set(df.cell)) == 10394
    assert set(df.bare.unique()) <= {0,1} and set(df.printed.unique()) <= {0,1}
    assert df.in_cluster.sum() == 1154
    summary = {}
    for name, known in [("published_followup_baseline",df.before.eq(1)),("all_verified",df.resolved.ge(0))]:
        assert not (known & df.resolved.lt(0)).any()
        summary[name], g = contrasts(df, known, df.source_stratum)
        if name == "all_verified":
            g.to_csv(HERE/"stratum_bounds.tsv", index=False, sep="\t")
    b = summary["published_followup_baseline"]
    assert b["n_verified"] == 1375 and b["partial"]["in_positive"] == 163 and b["partial"]["out_positive"]==487
    (HERE/"contrast_bounds.json").write_text(json.dumps(summary,indent=2)+"\n",encoding="utf8")
    print(json.dumps(summary["all_verified"],indent=2))


if __name__ == "__main__":
    main()
