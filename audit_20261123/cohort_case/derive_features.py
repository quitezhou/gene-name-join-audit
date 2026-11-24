"""Derive identity-discordant labels and diagnostic controls exclusively from feature tables."""
import csv
import gzip
import hashlib
from collections import defaultdict
from pathlib import Path

OUT = Path(__file__).resolve().parent
def feature_records():
    with (OUT/"deposit_inventory.tsv").open(encoding="utf8") as src:
        inventory = list(csv.DictReader(src, delimiter="\t"))
    panels = {}
    for r in inventory:
        if not r["path"].endswith("/features.tsv.gz"):
            continue
        data = gzip.decompress((OUT/"source"/(r["file_id"]+".tsv.gz")).read_bytes())
        fingerprint = hashlib.sha256(data).hexdigest()
        rows = [tuple(line.split("\t")) for line in data.decode().splitlines()]
        assert all(len(v)==3 and v[2]=="Gene Expression" for v in rows)
        panels[r["path"].rsplit("/",1)[0]] = (fingerprint, rows)
    return panels

def main():
    panels = feature_records()
    names = defaultdict(set)
    for rows in dict(panels.values()).values():
        for gene, symbol, _ in rows:
            names[symbol].add(gene)
    candidates = sorted(s for s, ids in names.items() if len(ids)>1)
    (OUT/"candidate_symbols.txt").write_text("\n".join(candidates)+"\n",encoding="utf8")
    controls = ["ACTB","GAPDH","B2M","EEF1A1","MALAT1","PTPRC","TMSB4X","RPLP0","RPL13A","UBC","TPT1","FTL","FTH1","HLA-A","RPS18","RPL10","RPS6","RPL7","RPL8","RPS4X","RPS3","RPS8","RPL11","RPL19","RPL21","RPS24","RPL27","RPL30","RPL31","RPS27","RPS28","RPS29","RPS9","RPS11","RPS14","RPS15A","RPS16","RPS19","RPS20","RPS23","RPL4","RPL5","RPL9","RPL12","RPL13","RPL14","RPL15","RPL17","RPL18","RPL22","RPL23","RPL24","RPL26","RPL29","RPL32","RPL34","RPL35","RPL36","RPL37","RPL38","RPL39","RPL41"]
    controls = [s for s in controls if len(names[s])==1 and all(any(row[1]==s for row in rows) for rows in dict(panels.values()).values())]
    (OUT/"control_symbols.txt").write_text("\n".join(controls)+"\n",encoding="utf8")

if __name__ == "__main__":
    main()

