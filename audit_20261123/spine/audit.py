"""Compare deposited indices; quantify spatial X on the nine discordant labels.

The external 10x matrix is used only as a label provenance control, after
checking its complete ordered ID vector against the spine portal exports.
"""
import argparse
import csv
import gzip
import hashlib
import json
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import h5py
import numpy as np
from scipy.sparse import csr_matrix
from fetch_axes import BASE, NAMES, OUT, RangeReader

HERE = OUT.parent
PROJ = HERE.parent.parent
def load_axis(name):
    meta = json.loads((OUT/(name+".json")).read_text())
    with gzip.open(OUT/(name+".tsv.gz"),"rt",encoding="utf-8",newline="") as f:
        rows = list(csv.DictReader(f,delimiter="\t"))
    ordered = [(r["stable_id"],r["index"]) for r in rows]
    digest = hashlib.sha256("".join(i+"\t"+s+"\n" for i,s in ordered).encode()).hexdigest()
    assert digest==meta["id_symbol_text_sha256"]
    return ordered, meta

def write_tsv(path, rows):
    with path.open("w",encoding="utf-8",newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t",lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

def inspect_matrices(lookup, meta):
    """Collect exported matrix values for discordant features, without refitting."""
    names = [r["spatial_label"] for r in lookup]
    idx = [r["spatial_col_0based"] for r in lookup]
    def scan(file):
        cache = OUT/(file+".expression.json")
        if cache.exists():
            cached=json.loads(cache.read_text())
            if cached["etag"]==meta[file]["etag"] and cached["labels"]==names:
                return cached["rows"]
        for attempt in range(3):
            try:
                with RangeReader(BASE+file+".h5ad") as stream, h5py.File(stream,"r") as f:
                    assert stream.etag==meta[file]["etag"]
                    # H5AD stores the matrix by row. Use larger requests only after opening it.
                    stream.cache.clear()
                    stream.size_chunk=8*1024**2
                    x=f["X"]
                    assert x.attrs["encoding-type"]=="csr_matrix"
                    mat=csr_matrix((x["data"][:],x["indices"][:],x["indptr"][:]),shape=tuple(x.attrs["shape"]))
                    select=mat[:,idx].tocsc()
                    select.sum_duplicates()
                    select.eliminate_zeros()
                    assert np.isfinite(select.data).all() and np.min(select.data,initial=0)>=0
                    signal=np.asarray(select.sum(axis=0,dtype=np.float64)).ravel()
                    rows=[{"spatial_file":file,"spatial_label":k,"n_spots":mat.shape[0],
                           "spots_nonzero":int(select.indptr[i+1]-select.indptr[i]),"sum_exported_X":float(signal[i])}
                          for i,k in enumerate(names)]
                cache.write_text(json.dumps({"etag":stream.etag,"labels":names,"rows":rows},indent=2)+"\n")
                print(file,"matrix summary complete",flush=True)
                return rows
            except Exception:
                if attempt==2: raise
    with ThreadPoolExecutor(3) as pool:
        matrices = list(pool.map(scan,NAMES[1:]))
    return [r for m in matrices for r in m]

def main(include_expression):
    sc_rows,sc_meta=load_axis(NAMES[0])
    sc=[s for i,s in sc_rows]
    assert all(not i for i,s in sc_rows) and len(sc)==len(set(sc))
    meta={}
    spatial=None
    for file in NAMES[1:]:
        rows,meta[file]=load_axis(file)
        if spatial is None:spatial=rows
        assert rows==spatial
    assert len(set(g for g,_ in spatial))==len(spatial) and len(set(s for _,s in spatial))==len(spatial)
    ctrl=PROJ/"cohort_direct_features/GSM6030901_patient_1_filtered_feature_bc_matrix.h5"
    with h5py.File(ctrl) as f:
        group=f["matrix/features"]
        comparator=list(zip(group["id"].asstr()[:],group["name"].asstr()[:]))
    assert [g for g,_ in comparator]==[g for g,_ in spatial]
    syms=Counter(s for _,s in comparator)
    sp_names={s for _,s in spatial}
    matched=set(sc)&sp_names
    discordant=[]
    for row,((gid,label),(_,base)) in enumerate(zip(spatial,comparator)):
        if label==base:continue
        assert syms[base]==2 and label==base+"-1"
        dot=base+".1"
        if dot not in sc:continue
        assert dot not in sp_names and label not in sc
        discordant.append({"base_symbol":base,"spatial_label":label,"atlas_label":dot,
                           "spatial_ensembl_id":gid,"spatial_col_0based":row,"atlas_row_0based":sc.index(dot),
                           "reference_id_status":"unavailable"})
    assert len(discordant)==len(sc)-len(matched)
    write_tsv(HERE/"discordant_labels.tsv",discordant)
    result={"pmcid":"PMC11579336","doi":"10.1038/s41467-024-54187-0","portal":"https://spinal-development.cellgeni.sanger.ac.uk/",
            "paper_call":"indeterminate","endpoint":"Direct exact-name overlap between public expression portal exports, before any hypothetical model fitting",
            "spatial_files":len(meta),"unique_ordered_spatial_id_name_vectors":1,"spatial_features":len(spatial),
            "spatial_spots_by_file":{k:v["n_obs"] for k,v in meta.items()},"atlas_features":len(sc),"atlas_cells":sc_meta["n_obs"],
            "direct_matches":len(matched),"discordant_atlas_labels":len(discordant),"candidate_overlap_after_targeted_aliases":len(matched)+len(discordant),
            "comparator":{"accession":"GSM6030901","url":"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6030901",
                          "local_file":ctrl.name,"sha256":hashlib.sha256(ctrl.read_bytes()).hexdigest(),
                          "ordered_stable_ids_identical":True,"different_final_names_at_spatial_positions":sum(a!=b for a,b in zip(spatial,comparator))},
            "interpretation":"The seven portal indices share the measured label discordance. The published model fit, retained features and signature IDs are unobserved in these objects."}
    (HERE/"result.json").write_text(json.dumps(result,indent=2)+"\n")
    if include_expression:
        rows=inspect_matrices(discordant,meta)
        write_tsv(HERE/"spatial_expression.tsv",rows)
        by_label = {}
        for label in [r['spatial_label'] for r in discordant]:
            part = [r for r in rows if r['spatial_label']==label]
            by_label[label] = {'spots_nonzero': sum(r['spots_nonzero'] for r in part),
                               'files_nonzero': sum(r['spots_nonzero']>0 for r in part),
                               'sum_exported_X': sum(r['sum_exported_X'] for r in part)}
        expression_summary = {'total_spots_across_files':sum(v['n_obs'] for v in meta.values()),'by_label':by_label}
        (HERE/'expression_summary.json').write_text(json.dumps(expression_summary,indent=2)+'\n')
        print('nonzero labels by file:',{k:sum(r["spots_nonzero"]>0 for r in rows if r["spatial_file"]==k) for k in meta},flush=True)

if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--expression', action='store_true', help='summarize deposited X for the discordant columns')
    main(parser.parse_args().expression)


