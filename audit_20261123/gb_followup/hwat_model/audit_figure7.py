"""Match both published Figure 7 output workbooks to native filtered spot barcodes."""
from concurrent.futures import ThreadPoolExecutor
import csv
import gzip
import sys
from pathlib import Path
import openpyxl
import numpy as np

ROOT=Path(__file__).resolve().parent
CASE=ROOT.parents[1]/"cohort_case"
sys.path.insert(0,str(CASE))
from cache_entry import cache

def main():
    inv={r["path"]:r for r in csv.DictReader((CASE/"deposit_inventory.tsv").open(),delimiter="\t")}
    groups=openpyxl.load_workbook(ROOT/"source/github_Figure_7a_and_ED_Figure_7a.xlsx",read_only=True,data_only=True)
    states=openpyxl.load_workbook(ROOT/"source/github_Figure_7cd_and_ED_Figure_7bc.xlsx",read_only=True,data_only=True)
    sheets=states.sheetnames
    def read_raw(s):
        key=f"{s}_trimmed/filtered_feature_bc_matrix/barcodes.tsv.gz"
        with gzip.open(cache(inv[key]),"rt") as src:
            return [b.strip() for b in src]
    with ThreadPoolExecutor(max_workers=5) as pool:
        raw=dict(zip(sheets,pool.map(read_raw,sheets)))
    rows=[]
    for s in sheets:
        gsheet=groups[s+"_cell2location"]
        fro=list(gsheet.values); to=list(states[s].values)
        gb=[r[0] for r in fro[1:]]; sb=[r[0] for r in to[1:]]; rb=raw[s]
        assert len(set(gb))==len(gb) and len(set(sb))==len(sb) and len(set(rb))==len(rb)
        vals=np.asarray([r[1:] for r in to[1:]],dtype=float)
        assert np.isfinite(vals).all() and vals.min()>=0
        rows.append([s,len(rb),len(gb),len(sb),len(fro[0])-1,len(to[0])-1,len(set(rb)&set(sb)),rb==sb,rb==gb,len(set(sb)&set(gb))])
    with (ROOT/"fig7_model_spot_axes.tsv").open("w",newline="") as out:
        w=csv.writer(out,delimiter="\t"); w.writerow(["section","raw_filtered_spots","group_output_spots","state_output_spots","group_columns","state_columns","raw_state_intersection","raw_state_order_identical","raw_group_order_identical","group_state_intersection"]);w.writerows(rows)
    print("Rows",sum(r[1] for r in rows),"intersections",sum(r[6] for r in rows))

if __name__=="__main__":
    main()
