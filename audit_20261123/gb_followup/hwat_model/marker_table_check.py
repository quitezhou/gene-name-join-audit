"""Trace preidentified names into the paper's supplementary marker lists.

For four released objects, export_marker_vectors.R preserves Idents and the
available RNA vectors. Compare their uncorrected detection fractions with the
reported statistics. A mismatch cannot be attributed to an identity error.
"""
import csv
from pathlib import Path
import pandas as pd
import openpyxl

from network_edges import FOCAL

HERE=Path(__file__).resolve().parent
OBJECTS={"T cells":"ly","B cells":"B","Endotelial":"endo","FAP_sc":"scFAP"}

def main():
    targets=FOCAL|{g+".1" for g in FOCAL}
    for num in [4,5]:
        w=openpyxl.load_workbook(HERE/f"source/MOESM{num}.xlsx",read_only=True,data_only=True)
        rows=[]
        for sheet in w:
            for i,r in enumerate(sheet.iter_rows(min_row=2,max_col=7,values_only=True),2):
                if r[0] in targets:
                    rows.append([sheet.title,i,*r])
        d=pd.DataFrame(rows,columns=["sheet","row","gene","p_val","avg_logFC","pct1","pct2","p_adj","cluster"])
        d.to_csv(HERE/f"MOESM{num}_hits.tsv",sep="\t",index=False)
    summaries=[]
    for sheet,obj in OBJECTS.items():
        f=pd.read_csv(HERE/f"source/{obj}_rows.tsv.gz",sep="\t")
        for r in d[d.sheet.str.strip().eq(sheet)].itertuples():
            pos=f.cluster==r.cluster
            known_cluster_count={"FAP_sc":17,"T cells":11,"B cells":2,"Endotelial":12}[sheet]
            same_count=known_cluster_count==f.cluster.nunique()
            pair={"object":obj,"sheet":sheet,"table_row":r.row,"gene":r.gene,"cluster":r.cluster,
                  "published_pct1":r.pct1,"published_pct2":r.pct2,"same_cluster_count":same_count}
            if pos.any() and same_count:
                v=f["counts_"+r.gene].gt(0)
                pair.update(rna_pct1=round(v[pos].mean(),3),rna_pct2=round(v[~pos].mean(),3))
                pair["both_percentages_match"]=pair["rna_pct1"]==r.pct1 and pair["rna_pct2"]==r.pct2
            summaries.append(pair)
    pd.DataFrame(summaries).to_csv(HERE/"marker_detection_validation.tsv",sep="\t",index=False)

if __name__=="__main__":
    main()
