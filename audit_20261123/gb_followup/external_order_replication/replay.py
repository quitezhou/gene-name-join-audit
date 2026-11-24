"""Replay the independent DLPFC comparison from exact section-level summaries."""
import argparse
import gzip
from itertools import combinations
import json
from pathlib import Path
import sys

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent / "systematic_mechanism"))
from audit_join import assign_keys
from verify_order import distinct_loci


def run(out):
    out.mkdir(parents=True,exist_ok=True)
    totals = pd.read_csv(HERE/"gene_totals.tsv.gz",sep="\t",keep_default_na=False).set_index("id")
    info = pd.read_csv(HERE/"sections.tsv",sep="\t",dtype={"section":str}).set_index("section")
    det = pd.read_csv(HERE/"duplicate_family_detections.tsv.gz",sep="\t",keep_default_na=False).set_index("id")
    with gzip.open(HERE.parents[2]/"gencode.v32.gene_metadata.json.gz","rt",encoding="utf8") as f:
        annotation = json.load(f)
    assert annotation["release"] == 32 and totals.index.is_unique
    full = dict(zip(assign_keys(totals.symbol, "-"),totals.index))
    assert len(full) == len(totals)
    complete = {id_:key for key,id_ in full.items()}
    axes = {}
    for s, row in info.iterrows():
        expressed = totals[totals[s] > 0]
        assert int(expressed[s].sum()) == row.total_umi
        axes[s] = dict(filtered=dict(zip(assign_keys(expressed.symbol,"-"),expressed.index)),
                       inherited={complete[i]:i for i in expressed.index},
                       by_id={i:i for i in expressed.index},unfiltered=full)
    pair_rows, records, strict_keys = [], [], set()
    for a,b in combinations(info.index,2):
        pair = dict(left=a,right=b,left_donor=info.at[a,"donor"],right_donor=info.at[b,"donor"])
        for arm in axes[a]:
            left,right = axes[a][arm],axes[b][arm]
            common = left.keys() & right.keys()
            count = strict = 0
            for key in sorted(common):
                first,second = left[key],right[key]
                if first == second:
                    continue
                marked = distinct_loci(annotation["gene_records"],first,second)
                records.append(dict(**pair,arm=arm,key=key,left_id=first,right_id=second,
                                    disjoint_v32=marked,left_umi=totals.at[first,a],right_umi=totals.at[second,b]))
                count += 1
                strict += int(marked)
                if marked:
                    strict_keys.add(key)
            pair_rows.append(dict(**pair,arm=arm,shared_keys=len(common),id_conflicts=count,disjoint=strict))
    pairs = pd.DataFrame(pair_rows)
    pairs.to_csv(out/"pairs.tsv",sep="\t",index=False)
    rec_cols = ["left","right","left_donor","right_donor","arm","key","left_id","right_id","disjoint_v32","left_umi","right_umi"]
    pd.DataFrame(records,columns=rec_cols).to_csv(out/"events.tsv",sep="\t",index=False)
    summary = dict(sections=len(info),donors=info.donor.nunique(),shared_complete_axis=len(totals),
                   spots=int(info.spots.sum()),umi=int(info.total_umi.sum()),strict_keys=sorted(strict_keys),arms={})
    for arm,group in pairs.groupby("arm"):
        within = group.left_donor == group.right_donor
        hit = group[group.disjoint > 0]
        summary["arms"][arm] = dict(all_pairs=len(group),id_positive=int((group.id_conflicts>0).sum()),
                                    disjoint_positive=len(hit),within_pairs=int(within.sum()),
                                    within_positive=int((within & (group.disjoint>0)).sum()),
                                    between_positive=int((~within & (group.disjoint>0)).sum()),
                                    sections_incident=len(set(hit.left)|set(hit.right)),
                                    donors_incident=len(set(hit.left_donor)|set(hit.right_donor)),
                                    largest_keys_in_pair=int(group.disjoint.max()))
    # A missing expressed-only key contributes zero. In the protected arm counts
    # always follow the complete-axis ID, including a genuine measured zero.
    layers = []
    for key in sorted(strict_keys):
        for s,observations in info.iterrows():
            for arm,id_ in [("filtered",axes[s]["filtered"].get(key)),("protected",full[key])]:
                row = dict(key=key,section=s,donor=observations.donor,arm=arm,id=id_ or "",
                           n_L3=observations.n_L3,n_WM=observations.n_WM)
                for label in ["L3","WM"]:
                    row["detected_"+label] = int(det.at[id_,s+"_"+label]) if id_ else 0
                row["contrast_pp"] = 100*(row["detected_L3"]/row["n_L3"] - row["detected_WM"]/row["n_WM"])
                layers.append(row)
    cols = ["key","section","donor","arm","id","n_L3","n_WM","detected_L3","detected_WM","contrast_pp"]
    layer = pd.DataFrame(layers,columns=cols)
    layer.to_csv(out/"layer_contrasts.tsv",sep="\t",index=False)
    summary["layer_donor_equal_weight"] = []
    if not layer.empty:
        donor = layer.groupby(["key","donor","arm"]).contrast_pp.mean().unstack()
        avg = donor.groupby("key").mean()
        avg["shift_pp"] = avg.filtered - avg.protected
        summary["layer_donor_equal_weight"] = avg.reset_index().to_dict("records")
        donor.reset_index().to_csv(out/"layer_donors.tsv",sep="\t",index=False)
    (out/"summary.json").write_text(json.dumps(summary,indent=2)+"\n",encoding="utf8")
    print(json.dumps(summary,indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--out",type=Path,default=HERE)
    run(parser.parse_args().out)
