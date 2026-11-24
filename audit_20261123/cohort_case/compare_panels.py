"""Compare symbol/ID dictionaries among submitted primary features, before any model filtering."""
from collections import defaultdict, Counter
import csv
from derive_features import feature_records, OUT
from trace_atlas import write_table

def main():
    panels=feature_records()
    summary=[]
    records=[]
    distinct={}
    for sample,(h,rows) in panels.items():
        names=defaultdict(list)
        for gene,symbol,_ in rows: names[symbol].append(gene)
        distinct[h]=names
        summary.append(dict(sample=sample,fingerprint=h,features=len(rows),unique_symbols=len(names),
                            duplicated_symbols=sum(len(v)>1 for v in names.values())))
    write_table("feature_fingerprints.tsv",summary)
    spat={h for sample,(h,_) in panels.items() if "trimmed" in sample}
    assert len(spat)==1
    spatial=distinct[next(iter(spat))]
    for h,ref in distinct.items():
        if h in spat:continue
        for symbol in sorted(spatial.keys()&ref.keys()):
            a,b=spatial[symbol],ref[symbol]
            if a[0]!=b[0]:
                records.append(dict(reference_fingerprint=h,symbol=symbol,spatial_ids="|".join(a),reference_ids="|".join(b),
                                    singleton_in_both=(len(a)==len(b)==1)))
    write_table("input_name_conflicts.tsv",records)

if __name__=="__main__":main()
