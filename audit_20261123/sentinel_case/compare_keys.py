"""Audit the seven published test rows, then compare names versus stable IDs.

Analysis operates on the deposited matrix rows. Stable-ID and pooled-symbol
arms only replace the seven preidentified collision keys in Supplement S11.
All other columns, QC, samples, libraries and weights are held fixed.
"""
import csv
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from anndata.utils import make_index_unique
from openpyxl import load_workbook

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]

def recorded_json(path):
    return json.loads(path.read_text(encoding='utf-8'))

def prepared():
    summaries=recorded_json(ROOT/'GSE278936_all_matrix_burden_result.json')
    sample_info={s['accession']:s for s in summaries['samples']}
    pub=pd.read_excel(HERE/'source/MOESM13.xlsx',index_col='gene').drop(columns='Unnamed: 0')
    targets=[k for k in summaries['affected_local_names'] if k in pub.index]
    pairs=recorded_json(ROOT/'prostate_spatial_versions_result.json')['comparison_after_per_table_anndata_uniquification']['mismatch_examples']
    pairs={p['local_unique_name']:p for p in pairs}
    info={s['accession']:s for s in recorded_json(HERE/'matrix_manifest.json')}
    s1=pd.read_excel(HERE/'source/MOESM3.xlsx',sheet_name='Discovery').set_index('Patient ID')
    archive=np.load(HERE/'pseudobulk_counts.npz')
    counts=[]
    bg=[]
    for acc,meta in sorted(info.items()):
        sample=sample_info[acc]
        label=sample['label']
        assert s1.loc[label.replace('_',' '),'Spots under tissue']==meta['raw_spots']==sample['spots']
        assert meta['raw_UMI']==sample['total_counts']
        path=ROOT/'panel_sources/all_features'/f'{acc}_{label}_features.tsv.gz'
        features=pd.read_csv(path,sep='\t',header=None,names=['id','symbol','type'])
        assert features['type'].eq('Gene Expression').all() and features.id.is_unique
        local=make_index_unique(pd.Index(features.symbol))
        symbol_to_row=pd.Series(np.arange(len(features)),index=local)
        id_to_row=pd.Series(np.arange(len(features)),index=features.id)
        old=sample['reference_cohort']=='GRCh38-3.0.0-content'
        qc=archive[f'{acc}_qc']
        raw=archive[f'{acc}_raw']
        assert len(qc)==len(features) and int(qc.sum())==meta['qc_UMI']
        bg.append(pd.Series(qc,index=local,name=label).reindex(pub.index, fill_value=0))
        for symbol,p in pairs.items():
            row=symbol_to_row[symbol]
            assert int(raw[row])==sample['counts_by_local_name'][symbol]
            group=features.symbol.eq(symbol).to_numpy()
            shared=set(p['first_symbol_group']) & set(p['second_symbol_group'])
            assert len(shared)==1
            anchor=next(iter(shared))
            stable_row=id_to_row[anchor]
            for stage,vector in [('raw',raw),('qc',qc)]:
                counts.append(dict(sample=label,group=label.split('_')[0],release=sample['reference_cohort'],symbol=symbol,
                            stage=stage,total_UMI=int(vector.sum()),local_id=features.id.iloc[row],stable_id=anchor,
                            local=int(vector[row]),stable=int(vector[stable_row]),pooled=int(vector[group].sum())))
    data=pd.DataFrame(counts)
    return pub,targets,data,pd.DataFrame(bg)

def main():
    pub,targets,data,bg=prepared()
    data.to_csv(HERE/'sample_key_counts.tsv',sep='\t',index=False)
    assert bg.shape==(39,17696) and data['sample'].nunique()==39
    out=[]
    for (stage,symbol),frame in data.groupby(['stage','symbol']):
        entry=dict(stage=stage,symbol=symbol)
        for arm in ('local','stable','pooled'):
            cpm=frame.assign(value=1e6*frame[arm]/frame.total_UMI).groupby('group').value.mean()
            entry[arm+'_TRNA_mean_CPM']=cpm['TRNA']
            entry[arm+'_NEADT_mean_CPM']=cpm['NEADT']
            entry[arm+'_log2_mean_CPM_ratio']=np.log2(cpm['NEADT']/cpm['TRNA']) if min(cpm)>0 else None
        out.append(entry)
    pd.DataFrame(out).to_csv(HERE/'paired_cpm_contrasts.tsv',sep='\t',index=False)
    pub.loc[targets].to_csv(HERE/'published_pseudobulk_rows.tsv',sep='\t')
    print(pd.DataFrame(out).query("stage=='qc'")[['symbol','local_log2_mean_CPM_ratio','stable_log2_mean_CPM_ratio','pooled_log2_mean_CPM_ratio']].to_string(index=False))
    # S4 -> S5 manual marker validation. Only the preidentified ten exact names
    # are included here, and no denominator is inferred from missing exports.
    s4=load_workbook(HERE/'source/MOESM6.xlsx',read_only=True,data_only=True)
    s5=load_workbook(HERE/'source/MOESM7.xlsx',read_only=True,data_only=True)
    summary=recorded_json(ROOT/'GSE278936_all_matrix_burden_result.json')
    samples={s['label'].replace('_',' '):s for s in summary['samples']}
    result=[]
    for region in s4.sheetnames:
        hits=[]
        for row in list(s4[region].values)[1:]:
            if row[1] in summary['affected_local_names']:
                sample=samples[row[-1]]
                hits.append((row[1],row[-1],sample['stable_id_by_affected_local_name'][row[1]]))
        assert len(hits)==len(set(hits))
        for row in list(s5[region].values)[1:]:
            symbol=row[1]
            if symbol in summary['affected_local_names'] and row[2]:
                those=[r for r in hits if r[0]==symbol]
                assert len(those)==row[2]
                result.append(dict(region=region,symbol=symbol,n_samples=row[2],stable_IDs=sorted({r[2] for r in those}),
                                   S4_sample_labels=[r[1] for r in those],S5_padj=row[3],S5_is_marker=row[4]=='x'))
    (HERE/'supplement_marker_trace.json').write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':
    main()
