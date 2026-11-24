"""Replay the published clinical signature and log-rank threshold series with both feature keys."""
import argparse
import hashlib
import json
from pathlib import Path
import zipfile

import numpy as np
import pandas as pd
from scipy.stats import chi2

from kidney_feature_audit import obtain, SOURCES

HERE = Path(__file__).resolve().parent


def logrank(time, event, low):
    observed_minus_expected = variance = 0.
    for t in np.unique(time[event==1]):
        cases = (time==t) & (event==1)
        risk = time>=t
        n, n1, d = risk.sum(), (risk&low).sum(), cases.sum()
        observed_minus_expected += (cases&low).sum() - d*n1/n
        if n > 1:
            variance += d*(n-d)*n1*(n-n1)/(n*n*(n-1))
    return float(chi2.sf(observed_minus_expected**2/variance,1))


def main(cache):
    SOURCES["clinical.zip"] = ("https://ndownloader.figshare.com/files/43043695","a8b1367a78b02fb673501de8d470170c")
    with zipfile.ZipFile(obtain(cache,"clinical.zip")) as z:
        with z.open("Source_data/Fig6/geneset_neptune.txt") as fh:
            bulk = pd.read_csv(fh,sep="\t")
        with z.open("Source_data/Fig6/clinical_data_neptune.csv") as fh:
            clinical = pd.read_csv(fh).set_index("SAFID")
    bulk = bulk.rename(columns={bulk.columns[0]:"stable_id",bulk.columns[1]:"symbol"})
    assert bulk.stable_id.is_unique and bulk.symbol.is_unique
    ref = pd.read_csv(HERE/"kidney_reference_axis.tsv.gz",sep="\t")
    sig = (HERE/"kidney_S10_symbols.txt").read_text(encoding="utf8").splitlines()
    trace=[]
    for s in sig:
        matches=ref[ref.symbol==s]
        if matches.empty:
            matches=ref[ref.symbol.str.upper()==s.upper()]
        assert len(matches)==1
        q=bulk[bulk.symbol==s]
        assert len(q)==1
        trace.append({"signature_name":s, "source_symbol":matches.symbol.iat[0],
                      "source_id":matches.stable_id.iat[0], "neptune_symbol":q.symbol.iat[0], "neptune_id":q.stable_id.iat[0]})
    tr=pd.DataFrame(trace)
    assert tr.source_id.is_unique
    tr.to_csv(HERE/"kidney_S10_ID_trace.tsv",sep="\t",index=False)
    patients=clinical.index.astype(str).tolist()
    assert clinical.index.is_unique and set(patients).issubset(bulk.columns)
    # Keep the original feature order in both arms to eliminate summation-order noise.
    nominal=bulk.loc[bulk.symbol.isin(sig),patients].mean(axis=0).to_numpy()
    shared=bulk.loc[bulk.stable_id.isin(tr.source_id),patients].mean(axis=0).to_numpy()
    scores=pd.DataFrame({"published_score":clinical.score.to_numpy(),"name_score":nominal,"id_score":shared})
    rank=np.argsort(scores.published_score.to_numpy())
    scores=scores.iloc[rank].reset_index(drop=True)
    scores.index.name="rank_zero_based"
    scores.to_csv(HERE/"kidney_Fig6_scores.tsv",sep="\t")
    s11=pd.read_csv(HERE/"kidney_S11_published.tsv",sep="\t")
    time=clinical.DaysBLtoESRDorEGFR40.to_numpy()
    event=clinical.ESRDorEGFR40.to_numpy()
    assert set(event)=={0,1} and (time>0).all()
    for name,values in [("name",nominal),("id",shared)]:
        thresholds=np.sort(values)[1:-1]
        assert len(np.unique(values)) == len(values) == len(s11)+2
        s11[name+"_threshold"]=thresholds
        s11[name+"_p"]=[logrank(time,event,values<t) for t in thresholds]
    s11.to_csv(HERE/"kidney_S11_replayed.tsv",sep="\t",index=False)
    result=dict(signature_genes=len(tr),direct_symbol_matches=int((tr.source_symbol==tr.signature_name).sum()),
                source_id_matches=int((tr.source_id==tr.neptune_id).sum()),clinical_scores=len(scores),
                maximum_score_residual=float(abs(scores.name_score-scores.published_score).max()),
                changed_scores=int((scores.name_score!=scores.id_score).sum()),thresholds=len(s11),
                maximum_threshold_residual=float(abs(s11.threshold-s11.name_threshold).max()),
                maximum_logrank_p_residual=float(abs(s11.p-s11.name_p).max()),
                changed_logrank_tests=int((s11.id_p!=s11.name_p).sum()),changed_p_lt_05=int(((s11.id_p<.05)!=(s11.p<.05)).sum()))
    print(result)
    assert result["maximum_score_residual"]<1e-12 and result["maximum_threshold_residual"]<1e-8
    assert result["maximum_logrank_p_residual"]<1e-7
    (HERE/"kidney_Fig6_summary.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf8")


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--cache",type=Path,required=True)
    main(p.parse_args().cache)
