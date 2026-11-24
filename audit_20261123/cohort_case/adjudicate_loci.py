"""Classify stable-ID pairs with the project's frozen same-assembly locus rule."""
import csv
import gzip
import json
from pathlib import Path
from derive_features import OUT
from trace_atlas import write_table

def read(n):
    with (OUT/n).open(encoding="utf8") as f:return list(csv.DictReader(f,delimiter="\t"))

def main():
    gtf=json.load(gzip.open(OUT.parent.parent/"gencode.v32.gene_metadata.json.gz","rt",encoding="utf8"))
    assert gtf["compressed_gtf_sha256"]=="0c36c37abb0c88637b04a8e401c93f0199a65bbdd0702c9e13bdfbadc053fe1e"
    genes=gtf["gene_records"]
    canonical={f"chr{i}" for i in range(1,23)}|{"chrX","chrY"}
    def locus(a,b):
        x,y=genes.get(a,[]),genes.get(b,[])
        if len(x)!=1 or len(y)!=1: return "unresolved"
        x,y=x[0],y[0]
        if {x["chrom"],y["chrom"]}=={"chrX","chrY"} or x["chrom"] not in canonical or y["chrom"] not in canonical:return "unresolved"
        if x["chrom"]!=y["chrom"] or x["end"]<y["start"] or y["end"]<x["start"]:return "disjoint"
        return "overlap"
    sp={r["output_label"]:r for r in read("spatial_S42_vectors.tsv")}
    data=[]
    for r in read("source_vectors.tsv"):
        s=sp.get(r["output_label"])
        if not s or s["stable_id"]==r["stable_id"]:continue
        data.append(dict(output_label=r["output_label"],atlas_object=r["object"],source=r["source"],atlas_id=r["stable_id"],
            spatial_id=s["stable_id"],atlas_raw_multiplicity=r["raw_name_multiplicity"],spatial_raw_multiplicity=s["raw_name_multiplicity"],
            atlas_match=r["match_mode"],atlas_cells=r["n_cells"],atlas_detected=r["detected_cells"],spatial_spots=s["n_spots"],
            spatial_detected=s["detected_spots"],locus_relation=locus(r["stable_id"],s["stable_id"])))
    write_table("cross_object_evidence.tsv",data)
    locations={}
    for r in read("confirmed_object_rows.tsv"):
        locations["within_atlas",r["output_label"],r["id_a"],r["id_b"]]=locus(r["id_a"],r["id_b"])
    for r in data:
        locations["across_deposits",r["output_label"],r["atlas_id"],r["spatial_id"]]=r["locus_relation"]
    result=[]
    for (kind,label,a,b),val in sorted(locations.items()):
        def fmt(i):
            z=genes.get(i,[])
            return "|".join(f'{q["chrom"]}:{q["start"]}-{q["end"]}({q["strand"]})' for q in z) or "absent"
        result.append(dict(comparison=kind,label=label,id_a=a,interval_a=fmt(a),id_b=b,interval_b=fmt(b),relation=val))
    write_table("locus_classification.tsv",result)
if __name__=="__main__":
    main()
