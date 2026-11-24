"""Verify archived Figure 3D and Supplementary Data 15 and the original vs repaired decisions."""
import hashlib
import json
from pathlib import Path
import tempfile
import numpy as np
import pandas as pd
import requests
from scipy.stats import chi2_contingency
from openpyxl import load_workbook
from fetch_axes import HERE, SRC

METHODS=['source','six_duplicate_groups_by_id','all_shared_ids']

def supplemental_sheet(i, sheet):
    name=f'41467_2024_48700_MOESM{i}_ESM.xlsx'
    manifest=json.loads((SRC/'supplement_workbooks.json').read_text())
    record=next(r for r in manifest if r['filename']==name)
    cache=Path(tempfile.gettempdir())/'nsclc_ref_check'/name
    if not cache.is_file():
        r=requests.get(record['url'],timeout=90);r.raise_for_status();cache.write_bytes(r.content)
    assert hashlib.sha256(cache.read_bytes()).hexdigest()==record['sha256']
    wb=load_workbook(cache,read_only=True,data_only=True)
    with_rows=wb[sheet].values; cols=next(with_rows); n=next((i for i,c in enumerate(cols) if c is None),len(cols))
    df=pd.DataFrame((r[:n] for r in with_rows if r[0] is not None),columns=cols[:n]);wb.close();return df

def main():
    fig=supplemental_sheet(6,'Fig3D').rename(columns={'Patient':'Patient_number','Label':'Tissue','LR':'pair'})
    e=pd.read_csv(HERE/'figure3d_endpoint.tsv',sep='\t');v=fig.merge(e,validate='one_to_one',on=['Patient_number','Tissue','pair']);assert len(v)==len(fig)==192
    v.to_csv(HERE/'figure3d_published_points.tsv',index=False,sep='\t');assert np.abs(v.source-v.coloc_pct).max()<1e-10
    tests=supplemental_sheet(5,'Supplementary Data 15').rename(columns={'Patient':'Patient_number','LR':'pair'});counts=pd.read_csv(HERE/'figure3d_spot_counts.tsv',sep='\t');data=[]
    for (pt,pair,method),paircounts in counts.groupby(['Patient_number','pair','method']):
        a=paircounts[['both','passing']].to_numpy();a[:,1]-=a[:,0]
        empty=np.any(a.sum(axis=0)==0);pvalue=np.nan if empty else chi2_contingency(a)[1]
        data.append(dict(Patient_number=pt,pair=pair,method=method,pvalue=pvalue,empty_margin=bool(empty)))
    data=pd.DataFrame(data).merge(tests,on=['Patient_number','pair'],validate='many_to_one');data['p_adj']=np.minimum(1,data.pvalue*(309*8));data['significant']=data.p_adj<.05
    base=data[data.method=='source']; assert len(base)==192 and base.empty_margin.sum()==2 and (base.loc[base.empty_margin,'padj']==1).all();assert (base.p_adj-base.padj).abs().max()<1e-10
    data.to_csv(HERE/'supp15_test_replay.tsv',index=False,sep='\t'); significant=data.pivot(index=['Patient_number','pair'],columns='method',values='significant')
    summary={'published_points':len(v),'maximum_baseline_residual_pp':float(np.abs(v.source-v.coloc_pct).max()),'published_tests':len(base),'calculable_tests':int((~base.empty_margin).sum()),'maximum_calculable_baseline_adjusted_p_residual':float((base.p_adj-base.padj).abs().max()),'empty_margin_published_p_is_one':int(base.empty_margin.sum()),'adjustment':309*8,'methods':{}}
    for method in METHODS:
        d=v[method]-v.source;n_pair=significant.groupby(level='pair')[method].sum(); summary['methods'][method]={'changed_points':int((abs(d)>1e-10).sum()),'largest_abs_change_pp':float(abs(d).max()),'changed_significance':int((significant[method]!=significant['source']).sum()),'significant_patient_by_pair_tests':int(significant[method].sum()),'pairs_with_four_or_more_patients':int((n_pair>=4).sum()),'pairs_passing_four_patients':n_pair[n_pair>=4].index.to_list()}
    (HERE/'figure3d_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf8');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
