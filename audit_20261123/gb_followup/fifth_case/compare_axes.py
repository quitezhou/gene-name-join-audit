"""Replay the published variable-axis operations for the named NSCLC files."""
import ast,csv,collections,gzip,json,re,warnings
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
from fetch_axes import SRC,HERE
def read_axis(key):
    a=pd.read_csv(SRC/(key+".tsv.gz"),sep="\t",names=["id","symbol","type"],dtype=str,keep_default_na=False)
    assert a.id.is_unique and a.type.eq("Gene Expression").all()
    return a

def key_map(a):
    var=ad.AnnData(shape=(0,len(a)))
    var.var_names=a.symbol.to_numpy()
    var.var_names_make_unique()
    b=a.copy();b.index=var.var_names
    return b

def loc(g1,g2,lookup):
    a,b=lookup.get(g1,[]),lookup.get(g2,[])
    if len(a)!=1 or len(b)!=1:return "unresolved"
    a,b=a[0],b[0];ch={a["chrom"],b["chrom"]}
    if not all(re.fullmatch(r"chr(?:[1-9]|1[0-9]|2[0-2]|X|Y)", c) for c in ch) or ch=={"chrX","chrY"}:return "unresolved"
    if len(ch)>1:return "disjoint"
    return "disjoint" if (a["end"]<b["start"] or b["end"]<a["start"]) else "overlap"

def main():
    config=ast.literal_eval(ast.parse((SRC/"config.py").read_text()).body[0].value)
    manifest=pd.read_csv(HERE/"input_axes.tsv",sep="\t",dtype={"accession_tail":str})
    raw=manifest[manifest.accession_tail=="526"]; vis=manifest[manifest.accession_tail=="530"]
    assert len(raw)==81 and len(vis)==40 and raw.axis_id.nunique()==vis.axis_id.nunique()==1
    r=read_axis(raw.axis_id.iat[0]);s=read_axis(vis.axis_id.iat[0])
    sk_id_ix=dict(zip(s.id,s.index))
    r_groups=r.groupby("symbol").id.apply(set);s_groups=s.groupby("symbol").id.apply(set)
    r=key_map(r)
    # Captioned C2L workflow concatenates the identical spatial axes by ID. The
    # 0.7.6 merge uses pandas intersect, which preserves this common order.
    s=key_map(s[~s.symbol.str.startswith("MT-")])
    refs=[pd.read_csv(f,sep="\t",keep_default_na=False).iloc[:,0].to_list() for f in sorted(SRC.glob("*.h5ad.var.tsv"))]
    assert len(refs)==2 and refs[0]==refs[1] and set(refs[0])<=set(r.index) and len(set(refs[0]))==len(refs[0])
    r=r.loc[refs[0]]; common=s.index[s.index.isin(r.index)]
    bad=[k for k in common if r.loc[k,"id"]!=s.loc[k,"id"]]
    gencode=json.load(gzip.open(HERE.parents[2]/"gencode.v32.gene_metadata.json.gz","rt"))["gene_records"]
    rows=[]
    for k in bad:
        ref_id,sp_id=r.loc[k,"id"],s.loc[k,"id"]
        rn,sn=r.loc[k,"symbol"],s.loc[k,"symbol"]
        rg,sg=r_groups[rn],s_groups[sn]
        kind="duplicate_associated" if max(len(rg),len(sg))>1 else "singleton_name"
        locus=loc(ref_id,sp_id,gencode)
        rows.append([k,rn,sn,ref_id,sp_id,len(rg),len(sg),bool(rg&sg),kind,locus])
    d=pd.DataFrame(rows,columns=["input_key","sc_base","vis_base","sc_source_id","vis_id","sc_base_members","vis_base_members","shared_id_in_source_groups","cause","gencode_v32_interval"])
    # Link each configured Sanger ID to the deposited sample identity by the SDRF.
    table=pd.read_csv(SRC/"E-MTAB-13530.sdrf.txt",sep="\t")
    mapping=table.assign(sanger=table["Scan Name"].str.extract(r"LNGsp(\d+)",expand=False)).groupby("sanger")["Source Name"].unique()
    assert all(len(v)==1 for v in mapping) and len(mapping)==40
    cohorts={}
    for env in ["Tumour","Background","Healthy"]:
        for value in config[env]["spatial"]:
            sample=mapping[re.search(r"LNGsp(\d+)",value)[1]][0]
            assert sample not in cohorts
            cohorts[sample]=env
    summary=[]
    for row in vis.itertuples(index=False):
        sample=row.file.split("-filtered")[0]
        z=np.load(SRC/(row.file+".npz"))
        for k,group in d.groupby(["cause","gencode_v32_interval"]):
            indices=[sk_id_ix[i] for i in group.vis_id]
            total=int(z["total"][indices].sum())
            summary.append((sample,cohorts[sample],row.spots,row.total_UMI,*k,len(indices),sum(z["total"][indices]>0),total))
    out=pd.DataFrame(summary,columns=["sample","environment","spots","total_UMI","cause","locus","different_ID_keys","expressed_keys","UMI_at_affected_columns"])
    out.to_csv(HERE/"per_slide.tsv",sep="\t",index=False)
    for k in ["Tumour","Background","Healthy"]:
        fit=[np.load(SRC/(row.file+".npz")) for row in vis.itertuples(index=False) if cohorts[row.file.split("-filtered")[0]]==k]
        d[k+"_UMI"]=[sum(int(z["total"][sk_id_ix[i]]) for z in fit) for i in d.vis_id]
        d[k+"_nonzero_slides"]=[sum(z["total"][sk_id_ix[i]]>0 for z in fit) for i in d.vis_id]
    d.to_csv(HERE/"discordant_model_keys.tsv",sep="\t",index=False)
    readouts=dict(collections.Counter(cohorts.values()))
    true=len(set(r.id)&set(s.id)); correct=len(common)-len(bad)
    out_summary=dict(reference_axis=len(r),nonMT_spatial_axis=len(s),matched_local=len(common),same_id=correct,different_id=len(bad),
                        shared_IDs=true,shared_identity_pairs_lost_by_local_join=true-correct,slides=readouts,
                        spots=int(vis.spots.sum()),UMI_all=int(vis.total_UMI.sum()),groups=[])
    lr_genes={g for pair in config["LR pairs"] for g in pair.split("_")}
    out_summary["configured_LR_pairs"]=len(config["LR pairs"])
    out_summary["configured_LR_gene_tokens"]=len(lr_genes)
    out_summary["different_ID_keys_in_configured_LR"]=sorted(lr_genes & set(bad))
    for (cause,locus),df in d.groupby(["cause","gencode_v32_interval"]):
        cols=["Tumour_UMI","Background_UMI","Healthy_UMI"]; counts=df[cols].sum(axis=1)
        out_summary["groups"].append(dict(cause=cause,locus=locus,matched_keys=len(df),nonzero=int((counts>0).sum()),UMI=int(counts.sum())))
    (HERE/"summary.json").write_text(json.dumps(out_summary,indent=2)+"\n",encoding="utf8")
    print(json.dumps(out_summary,indent=2))
if __name__=="__main__":main()
