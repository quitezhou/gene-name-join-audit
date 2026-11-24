"""Check how two exact-trace TSTD3 identifiers changed between fixed GENCODE releases."""
import csv
import gzip
import hashlib
import io
from pathlib import Path
import re

import requests

OUT=Path(__file__).resolve().parent
CHECKSUMS={27:"c1ed4f7a8b9bd868ef2f70ddfffa9360c627956e5fd35a2a5656d4e80e723c99",
           29:"c9bc82ee8b434089cc14ef4b8070dbea36654f8872c122f9ec3c253b010a6b3e",
           32:"0c36c37abb0c88637b04a8e401c93f0199a65bbdd0702c9e13bdfbadc053fe1e"}
IDS=("ENSG00000228439","ENSG00000279170")

def main():
    records=[]
    for v,expected_sha in CHECKSUMS.items():
        url=f"https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{v}/gencode.v{v}.primary_assembly.annotation.gtf.gz"
        response=requests.get(url,timeout=(15,120))
        response.raise_for_status()
        if hashlib.sha256(response.content).hexdigest()!=expected_sha:raise ValueError(url)
        found=set()
        with gzip.GzipFile(fileobj=io.BytesIO(response.content)) as data:
            for line in data:
                if b"\tgene\t" not in line or not any(i.encode() in line for i in IDS):continue
                chrom,_,_,start,end,_,strand,_,attrs=line.decode().rstrip().split("\t")
                attrs=dict(re.findall(r'(\w+) "([^"]+)"',attrs))
                stable=attrs["gene_id"].split(".")[0]
                if stable not in IDS:continue
                found.add(stable)
                records.append([v,stable,attrs["gene_id"],attrs["gene_name"],attrs["gene_type"],chrom,start,end,strand,url,expected_sha])
        for stable in sorted(set(IDS)-found):
            records.append([v,stable,"absent","","","","","","",url,expected_sha])
    with (OUT/"tstd3_history.tsv").open("w",encoding="utf8",newline="") as handle:
        writer=csv.writer(handle,delimiter="\t")
        writer.writerow("release stable_id versioned_id gene_name gene_type chrom start end strand gtf_url sha256".split())
        writer.writerows(records)
if __name__=="__main__":main()
