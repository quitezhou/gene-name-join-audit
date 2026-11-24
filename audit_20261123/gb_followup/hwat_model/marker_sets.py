"""Reconstruct Figure 1a marker Jaccard comparisons from Supplementary Data 1."""
from collections import defaultdict
from pathlib import Path
import csv
import gzip
import json

import numpy as np
import pandas as pd
import openpyxl
from scipy.sparse import csr_matrix

BASE = Path(__file__).resolve().parent
# Crosswalk from the authors' source-code abbreviations to the publication's sheet order.
COHORTS = ["J","K","P","H","N","O","G","A","F","D","Q","L","M","C","B","E","I"]

def read_marker_sets():
    records = []
    w = openpyxl.load_workbook(BASE/"source/MOESM4.xlsx", read_only=True, data_only=True)
    for label, s in zip(COHORTS, w):
        for row in s.iter_rows(min_row=2, max_col=7, values_only=True):
            gene,p,fc,up,down,adj,cluster = row
            if not isinstance(gene,str) or not isinstance(fc,(int,float)) or not isinstance(adj,(int,float)) or cluster is None:
                continue
            records.append((label, int(cluster),gene,fc,adj))
    return pd.DataFrame(records,columns=["cohort","cluster","gene","fc","adj"])

def similarities(frame):
    label = frame.cohort + "_" + frame.cluster.astype(str)
    ix, labs = pd.factorize(label)
    jx, features = pd.factorize(frame.gene)
    a = csr_matrix((np.ones(len(ix),dtype=np.int32),(ix,jx)),shape=(len(labs),len(features)))
    assert a.data.max()==1
    overlap = (a @ a.T).toarray()
    totals = np.asarray(a.sum(axis=1)).ravel()
    denom = totals[:,None]+totals[None,:]-overlap
    score = overlap/denom
    cohorts = np.asarray([lab.split("_")[0] for lab in labs])
    score[cohorts[:,None]==cohorts[None,:]]=np.nan
    return pd.DataFrame(score,index=labs,columns=labs)

def assess(d):
    orig = pd.read_csv(BASE/"source/Figure_1a_Jaccard.txt",sep=" ").groupby("cluster").jaccard.apply(list)
    out = []
    for key in orig.index.intersection(d.index):
        obs = sorted(orig.loc[key],reverse=True)[:5]
        calc = sorted(d.loc[key].dropna(),reverse=True)[:len(obs)]
        for rank,(published,reconstructed) in enumerate(zip(obs,calc),1):
            out.append(dict(cluster=key,rank=rank,published=published,reconstructed=reconstructed))
    out=pd.DataFrame(out)
    return out

def main():
    frame=read_marker_sets()
    frame.loc[frame.cohort.isin(["J","K","L","P","Q"])].to_csv(BASE/"source/eligible_marker_rows.tsv.gz",sep="\t",index=False)
    modes={}
    for mode,keep in {"up_adj":(frame.fc>0)&(frame.adj<.05),
                      "literal_code":((frame.fc>0)&(frame.adj<.05))|frame.cohort.isin(["J","K","L"])}.items():
        ranked=assess(similarities(frame.loc[keep]))
        ranked.to_csv(BASE/f"fig1_{mode}_audit.tsv",sep="\t",index=False)
        modes[mode] = dict(paired_rank_rows=len(ranked),exact=int(np.isclose(ranked.published,ranked.reconstructed,atol=1e-10,rtol=0).sum()),
                           median_absolute_difference=float(np.median(abs(ranked.published-ranked.reconstructed))))
    (BASE/"fig1_baseline_validation.json").write_text(json.dumps(modes,indent=2))

if __name__ == "__main__":
    main()
