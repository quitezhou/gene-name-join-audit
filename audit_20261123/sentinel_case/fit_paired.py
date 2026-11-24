"""Paired PyDESeq2 0.4.4 fits on the deposited-data reconstruction.

Only the seven S11 entries with different nominal IDs are substituted. The
17,689 background rows and their normalization factors are identical. The
other paper-specific annotation and pathology exclusions are unavailable.
"""
import argparse
import importlib.metadata as md
import json
from pathlib import Path
import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.default_inference import DefaultInference
from pydeseq2.preprocessing import deseq2_norm_fit,deseq2_norm_transform
from compare_keys import HERE,prepared

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('arm',choices=['local','stable','pooled'])
    p.add_argument('--cores',type=int,default=3)
    p.add_argument('--old-release',action='store_true',help='Post hoc within-release sensitivity on 17 TRNA and 11 NEADT specimens')
    args=p.parse_args()
    outdir=HERE/'old_release_sensitivity' if args.old_release else HERE
    outdir.mkdir(exist_ok=True)
    assert md.version('pydeseq2')=='0.4.4'
    pub,targets,data,baseline=prepared()
    assert baseline.shape==(39,17696)
    qc=data.query("stage=='qc'").pivot(index='sample',columns='symbol',values=args.arm)
    if args.old_release:
        old_samples=set(data.loc[data.release=='GRCh38-3.0.0-content','sample'])
        baseline=baseline.loc[baseline.index.isin(old_samples)]
    counts=baseline.copy()
    counts[targets]=qc.reindex(counts.index)[targets]
    assert baseline.drop(columns=targets).equals(counts.drop(columns=targets))
    if args.arm=='local':
        assert counts.equals(baseline)
    unaffected=baseline.drop(columns=targets).to_numpy()
    logmeans,keep=deseq2_norm_fit(unaffected)
    _,factors=deseq2_norm_transform(unaffected,logmeans,keep)
    metadata=pd.DataFrame({'treatment':[s.split('_')[0] for s in counts.index]},index=counts.index)
    assert metadata.treatment.value_counts().to_dict()=={'NEADT':11 if args.old_release else 22,'TRNA':17}
    infer=DefaultInference(n_cpus=args.cores)
    dds=DeseqDataSet(counts=counts,metadata=metadata,design_factors='treatment',ref_level=['treatment','TRNA'],inference=infer)
    dds.obsm['size_factors']=factors
    dds.layers['normed_counts']=dds.X / factors[:,None]
    # Follow 0.4.4 deseq2() verbatim after the common size-factor step.
    dds.fit_genewise_dispersions()
    dds.fit_dispersion_trend()
    dds.fit_dispersion_prior()
    dds.fit_MAP_dispersions()
    dds.fit_LFC()
    dds.calculate_cooks()
    if dds.refit_cooks:
        dds.refit()
    stats=DeseqStats(dds,inference=infer,contrast=['treatment','NEADT','TRNA'])
    stats.summary()
    out=stats.results_df
    out.to_csv(outdir/f'fit_{args.arm}.tsv.gz',sep='\t',compression='gzip')
    print(out.loc[targets].to_string())
    pd.DataFrame({'sample':counts.index,'size_factor':factors}).to_csv(outdir/f'size_factors_{args.arm}.tsv',sep='\t',index=False)
    versions={k:md.version(k) for k in ['pydeseq2','numpy','pandas','scipy','anndata','statsmodels','scikit-learn']}
    (outdir/f'fit_versions_{args.arm}.json').write_text(json.dumps(versions,indent=2)+'\n')

if __name__=='__main__':
    main()

