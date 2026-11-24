"""Check the derived human-kidney feature and clinical results without access to a data server."""
import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent


def tab(name):
    return pd.read_csv(HERE/name, sep="\t")


def main():
    ref = tab("kidney_reference_axis.tsv.gz")
    vis = tab("kidney_spatial_axis.tsv.gz")
    removed = ref[ref.symbol.duplicated()]
    assert len(ref)==23446 and len(vis)==17943 and ref.stable_id.is_unique and vis.stable_id.is_unique
    assert len(removed)==116 and removed.symbol.nunique()==42 and removed.umi.sum()==163532
    assert (removed.cells>=50).all() and (removed.cells>=50+2623).sum()==7
    audited = vis.drop_duplicates("symbol").merge(ref.drop_duplicates("symbol"),on="symbol",suffixes=("_sp","_ref"))
    assert len(audited)==14336
    mis=audited[audited.stable_id_sp!=audited.stable_id_ref].set_index("symbol").sort_index()
    recorded=tab("kidney_shared_name_mismatches.tsv").set_index("symbol").sort_index()
    assert len(mis)==11 and mis.index.equals(recorded.index)
    for left,right in [("stable_id_sp","spatial_id"),("stable_id_ref","ref_first_id"),("umi1","spatial1_umi"),
                       ("umi2","spatial2_umi"),("umi","reference_first_umi"),("cells","reference_first_detected")]:
        assert np.array_equal(mis[left],recorded[right])
    groups=recorded[recorded.ref_symbol_rows>1]
    assert len(groups)==5 and (groups.spatial_symbol_rows==1).all()
    dropped=removed[removed.stable_id.isin(groups.spatial_id)]
    assert set(dropped.loc[dropped.cells>=50+2623,"symbol"])=={"ATXN7","SCO2"}
    tr=tab("kidney_S10_ID_trace.tsv")
    syms=(HERE/"kidney_S10_symbols.txt").read_text().splitlines()
    assert tr.signature_name.tolist()==syms and tr.signature_name.is_unique and len(tr)==65
    assert (tr.source_id==tr.neptune_id).all() and tr.source_id.is_unique
    assert (tr.source_symbol==tr.signature_name).sum()==64
    assert (tr.signature_name.str.upper()==tr.source_symbol.str.upper()).all()
    for row in tr.itertuples():
        candidates=ref[ref.symbol.str.upper()==row.signature_name.upper()]
        assert len(candidates)==1 and candidates.stable_id.iat[0]==row.source_id
    scores=tab("kidney_Fig6_scores.tsv")
    assert len(scores)==134 and (scores.name_score==scores.id_score).all()
    assert abs(scores.name_score-scores.published_score).max()<1e-12
    tests=tab("kidney_S11_replayed.tsv")
    pub=tab("kidney_S11_published.tsv")
    assert np.allclose(tests.p,pub.p) and np.array_equal(tests.n_low,np.arange(1,133))
    assert np.array_equal(tests.n_high,134-tests.n_low)
    assert np.allclose(tests.threshold,np.sort(scores.name_score)[1:-1],rtol=0,atol=1e-8)
    assert np.array_equal(tests.name_p,tests.id_p) and abs(tests.p-tests.name_p).max()<1e-7
    print("Kidney source features and clinical score replay verified")


if __name__=="__main__":
    main()
