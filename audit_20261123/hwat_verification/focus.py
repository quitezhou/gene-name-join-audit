"""Pick the source and object witnesses referred to in the independent report."""
from verify import HERE,CASE,tab,write
def main():
    all_rows=tab(HERE/"atlas_all_rows.tsv")
    all_spatial=tab(HERE/"spatial_all_rows.tsv")
    within=tab(CASE/"locus_classification.tsv")
    switched={r["label"] for r in within if r["comparison"]=="within_atlas"}
    out=[]
    for r in all_rows:
        if r["anchor"]!="barcode_and_metadata" or r["label"] not in switched or r["object"]!="adipo":continue
        if r["source"] not in {"Massier et al. (#1)","Massier et al. (#4)/21L005051"}:continue
        assert r["expected_id"]==r["support_id"] and r["n_support_rows"]=="1"
        out.append(r)
    assert len(out)==10
    write("within_adipo_verified.tsv",out)
    sp={r["label"]:r for r in all_spatial}
    cross=[r for r in within if r["comparison"]=="across_deposits" and r["relation"]=="disjoint"]
    out=[]
    for pair in cross:
        label=pair["label"]
        obj="endo" if label=="AP001528.1" else "ly"
        l=next(r for r in all_rows if r["source"]=="Massier et al. (#4)/21L005051" and r["anchor"]=="barcode_and_metadata" and r["object"]==obj and r["label"]==label)
        s=sp[label]
        assert l["support_id"]==pair["id_a"] and s["exact_id"]==pair["id_b"] and l["n_support_rows"]=="1" and s["n_exact_rows"]=="1"
        out.append({"label":label,"atlas_object":obj,"atlas_source":l["source"],"atlas_n_cells":l["matched_barcodes"],
                    "atlas_n_positive":l["target_n_detected"],"atlas_id":l["support_id"],"atlas_competing_raw_features":l["raw_features"],
                    "atlas_best_other_support_jaccard":l["best_other_support_jaccard"],"spatial_sample":"S42","spatial_id":s["exact_id"],
                    "spatial_n_spots":s["matched_barcodes"],"spatial_n_positive":s["target_n_detected"],"spatial_competing_raw_features":s["raw_features"],
                    "interval_a":pair["interval_a"],"interval_b":pair["interval_b"],"relation":pair["relation"]})
    assert len(out)==3
    write("cross_disjoint_verified.tsv",out)

if __name__=="__main__":main()
