"""Compare detection fractions for a bare symbol versus a traced shared ID.

Only nonzero rows uniquely traced by the directly extracted adipocyte dataset
enter this comparison. Compute source-specific fractions on the same cells in
both interpretations. This is a descriptive reuse example, not a paper endpoint.
"""
from verify import HERE,tab,write

def main():
    rows=tab(HERE/"base_and_suffix_trace.tsv")
    sources=["Massier et al. (#1)","Massier et al. (#4)/21L005051"]
    results=[]
    for symbol,shared in [("PINX1","ENSG00000254093"),("ARMCX5-GPRASP2","ENSG00000286237")]:
        fractions={}
        summary={"symbol":symbol,"shared_id":shared}
        for i,s in enumerate(sources,1):
            family=[r for r in rows if r["sample"]==s and r["stored_k"] in [symbol,symbol+".1"]]
            bare=next(r for r in family if r["stored_k"]==symbol)
            match=[r for r in family if r["unique_id"]==shared]
            assert len(match)==1 and bare["support_matches"]=="1" and match[0]["support_matches"]=="1"
            traced=match[0]
            n=int(bare["n_cells"])
            assert n==int(traced["n_cells"])
            for kind,x in [("bare",bare),("shared_id",traced)]:
                fractions[kind,i]=100*int(x["n_detected"])/n
                summary[f"{kind}_positive_{i}"]=int(x["n_detected"])
                summary[f"{kind}_pct_{i}"]=round(fractions[kind,i],6)
            summary[f"cells_{i}"]=n
            summary[f"shared_id_row_{i}"]=traced["stored_k"]
        for kind in ("bare","shared_id"):
            summary[kind+"_pct_point_change"]=round(fractions[kind,2]-fractions[kind,1],6)
        results.append(summary)
    write("adipocyte_detection_contrasts.tsv",results)

if __name__=="__main__":main()
