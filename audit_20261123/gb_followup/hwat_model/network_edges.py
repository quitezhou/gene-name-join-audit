"""Compare archived Figure 1b edge weights with the published marker tables.

Keep cluster markers fixed. Separately recalculate overlaps for the five traced
families, for all families duplicated in either panel, and for unambiguous
singleton name aliases. The last overlay counts same-name/different-ID losses
and missing pairs that have identical IDs under two distinct unused names.
"""
from collections import Counter, defaultdict
import csv
import io
from pathlib import Path
import sys
import zipfile

import pandas as pd

HERE=Path(__file__).resolve().parent
STUDY=HERE.parents[1]/"cohort_case"
sys.path.insert(0,str(STUDY))
from derive_features import feature_records

FOCAL={"AHRR","ARMCX5-GPRASP2","ATXN7","PINX1","POLR2J3"}

def cytab(z,kind):
    name=next(n for n in z.namelist() if "/18153-" in n and "SHARED_ATTRS-" in n and "Cy"+kind in n)
    rows=list(csv.reader(io.StringIO(z.read(name).decode("utf-8"))))
    return pd.DataFrame(rows[5:],columns=rows[1])

def name_table(rows):
    # R-style name repair, respecting existing literal dotted symbols.
    reserve={s for i,s,t in rows}
    seen=defaultdict(int)
    out={}
    for gene,s,_ in rows:
        if not seen[s]:
            label=s
        else:
            while f"{s}.{seen[s]}" in reserve:
                seen[s]+=1
            label=f"{s}.{seen[s]}"
        assert label not in out
        out[label]=(gene,s)
        seen[s]+=1
    return out,Counter(s for _,s,_ in rows)

def main():
    all_panels=feature_records()
    # Authors' Jaccard code maps J=lipidlab (#1), K=leipzig_sc (#2),
    # L=leipzig_om (#4), P=lipid_1s (#3), Q=Backdahl spatial.
    cohorts={"J":"Massier et al. (#1)", "K":"Massier et al. (#2)/21L005039",
             "L":"Massier et al. (#4)/21L005051","P":"Massier et al. (#3)",
             "Q":next(k for k in all_panels if k.startswith("S42"))}
    panels={k:name_table(all_panels[v][1]) for k,v in cohorts.items()}
    f=pd.read_csv(HERE/"source/eligible_marker_rows.tsv.gz",sep="\t")
    f=f[((f.fc>0)&(f.adj<.05))|f.cohort.isin(["J","K","L"])]
    sets={f"{c}_{i}":set(g.gene) for (c,i),g in f.groupby(["cohort","cluster"])}
    with zipfile.ZipFile(HERE/"source/Figure_1b.cys") as z:
        nodes=cytab(z,"Node").set_index("id")
        edges=cytab(z,"Edge")
    output=[]
    for _,e in edges.iterrows():
        a,b=nodes.loc[e.source,"label"],nodes.loc[e.target,"label"]
        if a not in sets or b not in sets:
            continue
        left,right=sets[a],sets[b]
        ma,na=panels[a[0]]; mb,nb=panels[b[0]]
        ol=left&right
        loss=sorted(s for s in ol if s in ma and s in mb and na[s]==nb[s]==1 and ma[s][0]!=mb[s][0])
        span={s for s,n in na.items() if n>1}|{s for s,n in nb.items() if n>1}
        def targeted_key(s,m,families):
            if s in m and m[s][1] in families:
                return "verified_family_id:"+m[s][0]
            return s
        lnew={targeted_key(s,ma,FOCAL) for s in left}
        rnew={targeted_key(s,mb,FOCAL) for s in right}
        ldup={targeted_key(s,ma,span) for s in left}
        rdup={targeted_key(s,mb,span) for s in right}
        assert len(ldup)==len(left) and len(rdup)==len(right)
        assert len(lnew)==len(left) and len(rnew)==len(right)
        la={ma[s][0]:s for s in left if s in ma and na[ma[s][1]]==1}
        rb={mb[s][0]:s for s in right if s in mb and nb[mb[s][1]]==1}
        gained=sorted((gene,la[gene],rb[gene]) for gene in la.keys() & rb.keys()
                      if la[gene]!=rb[gene] and la[gene] not in right and rb[gene] not in left)
        weight=int(e.weight)
        actual=len(ol)
        eq=actual==weight
        # Verify the individual exported denominators through both percentages.
        left_exact=abs(float(e.Percent_1)-100*actual/len(left))<1e-3
        right_exact=abs(float(e.Percent_2)-100*actual/len(right))<1e-3
        output.append(dict(left=a,right=b,source_cell=nodes.loc[e.source,"annotation"],target_cell=nodes.loc[e.target,"annotation"],
                           paper_overlap=weight,reconstructed_overlap=actual,left_size=len(left),right_size=len(right),
                           weight_and_both_percentages_match=eq and left_exact and right_exact,
                           focal_overlap=len(lnew&rnew),focal_delta=len(lnew&rnew)-actual,
                           all_duplicate_family_overlap=len(ldup&rdup),all_duplicate_family_delta=len(ldup&rdup)-actual,
                           identical_singleton_label_different_id=len(loss),singleton_labels="|".join(loss),
                           recovered_singleton_alias_pairs=len(gained),
                           unambiguous_singleton_overlay=actual-len(loss)+len(gained),
                           alias_traces="|".join(",".join(x) for x in gained)))
    d=pd.DataFrame(output)
    d.to_csv(HERE/"fig1b_tested_edges.tsv",sep="\t",index=False)
    ok=d[d.weight_and_both_percentages_match]
    print(len(edges),len(d),len(ok),'focal',ok.focal_delta.value_counts().to_dict(),
          'singleton',ok.identical_singleton_label_different_id.value_counts().to_dict())
    print(ok[ok.focal_delta!=0].head(10).to_string(index=False))
    print("alias recoveries",ok.recovered_singleton_alias_pairs.describe())

if __name__=="__main__":
    main()
