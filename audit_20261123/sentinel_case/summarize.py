"""Validate all three fits, summarize them and render a figure from the outputs."""
import csv
import gzip
import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from compare_keys import HERE,ROOT,prepared,recorded_json

def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    pubs,genes,_,bg=prepared()
    dfs={'published':pubs}
    qc=recorded_json(HERE/'matrix_manifest.json')
    assert len(qc)==39 and len({m['accession'] for m in qc})==39
    assert not recorded_json(HERE/'unavailable.json')
    assert len({digest(HERE/f'size_factors_{arm}.tsv') for arm in ['local','stable','pooled']})==1
    for arm in ['local','stable','pooled']:
        df=pd.read_csv(HERE/f'fit_{arm}.tsv.gz',sep='\t',index_col=0)
        assert df.index.equals(pubs.index)
        dfs[arm]=df
    allstats=pd.concat({k:v.loc[genes,['log2FoldChange','lfcSE','padj']] for k,v in dfs.items()},axis=1)
    allstats.columns=['_'.join(c) for c in allstats.columns]
    allstats.to_csv(HERE/'paired_DE_results.tsv',sep='\t')
    loci=json.load(gzip.open(ROOT/'gencode.v32.gene_metadata.json.gz'))['gene_records']
    pairs=recorded_json(ROOT/'prostate_spatial_versions_result.json')['comparison_after_per_table_anndata_uniquification']['mismatch_examples']
    mappings=[]
    for pair in pairs:
        symbol=pair['local_unique_name']
        id1,id2=pair['first_stable_id'],pair['second_stable_id']
        def loc(ID):
            hits=loci.get(ID,[])
            return hits[0] if len(hits)==1 else None
        l1,l2=loc(id1),loc(id2)
        def readable(l):
            return f"{l['chrom']}:{l['start']}-{l['end']}({l['strand']})" if l else 'unresolved'
        cls='unresolved'
        if l1 and l2:
            overlap=(l1['chrom']==l2['chrom'] and min(l1['end'],l2['end'])>=max(l1['start'],l2['start']))
            cls=('overlap_same_strand' if l1['strand']==l2['strand'] else 'overlap_opposite_strand') if overlap else 'nonoverlap'
        first=set(pair['first_symbol_group'])
        second=set(pair['second_symbol_group'])
        mappings.append(dict(gene=symbol,in_S11=symbol in genes,old_unsuffixed_id=id1,new_unsuffixed_id=id2,
                            old_group=','.join(sorted(first)),new_group=','.join(sorted(second)),
                            common_id=next(iter(first&second)),old_locus_v32=readable(l1),new_locus_v32=readable(l2),
                            coordinate_class=cls))
    pd.DataFrame(mappings).to_csv(HERE/'ID_mapping.tsv',sep='\t',index=False)
    sub=pd.concat([dfs[k].log2FoldChange.rename(k) for k in ['published','local']],axis=1).replace([np.inf,-np.inf],np.nan).dropna()
    def hits(df):
        return set(df.index[(df['padj']<.05) & (df.log2FoldChange.abs()>=1)])
    called={a:hits(df) for a,df in dfs.items()}
    summary=dict(samples=len(qc),raw_spots=sum(m['raw_spots'] for m in qc),qc_spots=sum(m['qc_spots'] for m in qc),
                 tested_panel=len(bg.columns),unchanged_background=len(bg.columns)-len(genes),studied_entries=genes,
                 finite_paired_baseline_LFC=len(sub),all_genes_baseline_LFC_spearman=spearmanr(sub.published,sub.local).statistic,
                 all_genes_baseline_LFC_median_abs_difference=(sub.published-sub.local).abs().median(),
                 all_genes_baseline_LFC_90pct_abs_difference=(sub.published-sub.local).abs().quantile(.9),
                 candidates_baseline_LFC_max_abs_difference=(sub.loc[genes,'published']-sub.loc[genes,'local']).abs().max(),
                 paper_cutoff='abs(log2FC)>=1 and adjusted p<0.05',all_genes_paper_hits=len(called['published']),
                 all_genes_local_hits=len(called['local']),all_genes_published_local_hits_intersection=len(called['published']&called['local']),
                 candidate_hits={a:sorted(set(genes)&h) for a,h in called.items()})
    (HERE/'validation.json').write_text(json.dumps(summary,indent=2)+'\n')
    import ast
    nb=recorded_json(HERE/'source'/'spatial-pseudobulk-data-analysis.ipynb')
    def gene_list(cell):
        tree=ast.parse(''.join(nb['cells'][cell]['source']))
        group=next(x for x in tree.body if isinstance(x,ast.Assign) and any(isinstance(y,ast.Name) and y.id=='genes' for y in x.targets))
        return [s for s in ast.literal_eval(group.value) if s.strip()]
    g5f,g5i=gene_list(19),gene_list(21)
    assert not (set(g5f+g5i) & set(genes))
    def stars(p):
        return '***' if p<.001 else '**' if p<.01 else '*' if p<.05 else 'none'
    fig5i=pd.concat({k:v.loc[g5i,'padj'] for k,v in dfs.items()},axis=1)
    for arm in dfs:
        fig5i[arm+'_stars']=fig5i[arm].map(stars)
    fig5i.to_csv(HERE/'fig5i_marker_check.tsv',sep='\t')
    summary.update(fig5f_listed_genes=len(g5f),fig5i_listed_genes=len(g5i),
                   changed_fig5i_stars_local_vs_stable=int((fig5i.local_stars!=fig5i.stable_stars).sum()),
                   changed_fig5i_stars_local_vs_published=int((fig5i.local_stars!=fig5i.published_stars).sum()))
    (HERE/'validation.json').write_text(json.dumps(summary,indent=2)+'\n')
    sources=HERE/'source'
    records={}
    for p in sources.iterdir():
        if p.suffix=='.xlsx':
            idx=p.stem.replace('MOESM','')
            url=f'https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-024-54364-1/MediaObjects/41467_2024_54364_MOESM{idx}_ESM.xlsx'
        else:
            url='https://raw.githubusercontent.com/akiviaho/ST-prostate/df0d2261af3a14a3fa0c65b553a200d3edcbe279/'+p.name
        records[p.name]={'url':url,'bytes':p.stat().st_size,'sha256':digest(p)}
    (HERE/'primary_sources.json').write_text(json.dumps(records,indent=2)+'\n')
    fig,ax=plt.subplots(figsize=(8.1,4.9),layout='constrained')
    y=np.arange(len(genes))[::-1]
    for v,z in zip(genes,y):
        row=allstats.loc[v]
        ax.plot([row.local_log2FoldChange,row.stable_log2FoldChange],[z,z],color='#9ba8b2',lw=1.8,zorder=1)
    ax.axvline(0,color='#6e7781',lw=.8)
    ax.scatter(dfs['local'].loc[genes,'log2FoldChange'],y,marker='o',facecolors='white',edgecolors='#222e3a',s=54,label='Local symbol',zorder=3)
    ax.scatter(dfs['stable'].loc[genes,'log2FoldChange'],y,marker='o',color='#04765e',s=28,label='Shared stable ID',zorder=4)
    ax.scatter(dfs['published'].loc[genes,'log2FoldChange'],y+.15,marker='x',color='#54627d',s=31,label='Published S11',zorder=2)
    ax.set_yticks(y,genes);ax.set_xlim(-3.4,11.45);ax.set_ylim(-.6,6.65)
    ax.set_xlabel('Treatment contrast (NEADT versus TRNA), log₂ fold change')
    ax.set_title('Gene identity and the prostate pseudobulk contrast',loc='left',weight='bold',pad=12)
    ax.legend(loc='upper center',bbox_to_anchor=(.5,-.17),frameon=False,ncol=3)
    ax.spines[['top','right','left']].set_visible(False)
    ax.tick_params(axis='y',length=0)
    ax.xaxis.grid(color='#e5e8eb',lw=.6)
    for ext in ['png','pdf']:
        fig.savefig(HERE/f'prostate_pseudobulk_ID_repair.{ext}',dpi=200,bbox_inches='tight')
    plt.close(fig)

if __name__=='__main__':
    main()

