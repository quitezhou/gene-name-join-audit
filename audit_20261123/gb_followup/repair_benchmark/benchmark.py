"""Compare observed name-based joins with literal source-ID pairs on deposited axes."""
import argparse
import csv
import gzip
import hashlib
import itertools
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PANEL = HERE.parent / "systematic_mechanism"
CASE = HERE.parent / "fifth_case"
sys.path.insert(0, str(PANEL))
from audit_join import Axis, load
from check_lung_join import records, source_axis


def write_table(out, name, rows):
    with (out / name).open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def repeated(rows):
    return {symbol for symbol, n in Counter(r.symbol for r in rows).items() if n > 1}


def analyze(left, right, source_duplicates, amount=None, inherited=None):
    # All supplied axes must have unique source IDs and keys; ambiguous IDs are not silently collapsed.
    left_by_id, right_by_id = ({r.id: r for r in side} for side in (left, right))
    left_keys, right_keys = ({r.key: r for r in side} for side in (left, right))
    assert len(left) == len(left_by_id) == len(left_keys)
    assert len(right) == len(right_by_id) == len(right_keys)
    shared = left_by_id.keys() & right_by_id.keys()
    local_pairs = [(left_keys[k], right_keys[k]) for k in left_keys.keys() & right_keys.keys()]
    final_duplicates = repeated(left) | repeated(right)
    def drop(group):
        return [(l, r) for l, r in local_pairs if l.symbol not in group and r.symbol not in group]
    arms = {"local_keys": local_pairs, "drop_source_duplicates": drop(source_duplicates),
            "drop_final_duplicates": drop(final_duplicates),
            "exact_id": [(left_by_id[i], right_by_id[i]) for i in shared]}
    if inherited is not None:
        # The caller has established ordered subset identity; inherit the master's original keys.
        master_left, master_right = ({inherited[r.id]: r for r in side} for side in (left, right))
        arms["inherited_master_keys"] = [(master_left[k], master_right[k]) for k in master_left.keys() & master_right.keys()]
    local_correct = {l.id for l, r in local_pairs if l.id == r.id}
    summary = []
    for arm, pairs in arms.items():
        correct = {l.id for l, r in pairs if l.id == r.id}
        wrong = [(l, r) for l, r in pairs if l.id != r.id]
        wrong_spatial = {r.id for l, r in wrong}
        dup_wrong = {r.id for l, r in wrong if l.symbol in source_duplicates or r.symbol in source_duplicates}
        assert len(pairs) == len(correct) + len(wrong_spatial) and correct <= shared
        row = dict(arm=arm, left_rows=len(left), right_rows=len(right), exact_shared_ids=len(shared),
                   left_only_ids=len(left)-len(shared), right_only_ids=len(right)-len(shared),
                   matched_pairs=len(pairs), same_id_pairs=len(correct), different_id_pairs=len(wrong),
                   wrong_source_duplicate=len(dup_wrong), wrong_source_singleton=len(wrong)-len(dup_wrong),
                   shared_ids_missed=len(shared-correct), correct_coverage_pct=round(100*len(correct)/len(shared), 5),
                   correct_ids_lost_vs_local=len(local_correct-correct), correct_ids_gained_vs_local=len(correct-local_correct))
        if amount is not None:
            mol = lambda ids: sum(amount[i] for i in ids)
            row.update(spatial_umi_analyzed=mol(right_by_id), spatial_umi_shared_ids=mol(shared),
                       spatial_umi_right_only_ids=mol(right_by_id.keys()-shared), spatial_umi_correct=mol(correct),
                       spatial_umi_at_wrong_pairs=mol(wrong_spatial), spatial_umi_at_wrong_source_duplicate=mol(dup_wrong),
                       spatial_umi_at_wrong_source_singleton=mol(wrong_spatial-dup_wrong),
                       spatial_umi_on_missed_shared_ids=mol(shared-correct), spatial_umi_correct_lost_vs_local=mol(local_correct-correct),
                       spatial_umi_correct_gained_vs_local=mol(correct-local_correct))
        summary.append(row)
    return summary


def lung(out):
    manifest = records(CASE / "input_axes.tsv")
    sc, _ = source_axis(manifest, "526", 81)
    spatial, _ = source_axis(manifest, "530", 40)
    full_sc, full_spatial = (list(Axis(entries).by_key.values()) for entries in (sc, spatial))
    assert all(len({r.id for r in axis}) == len(axis) for axis in (full_sc, full_spatial))
    full_sc_by_key = {r.key: r for r in full_sc}
    var = [[r["_index"] for r in records(f)] for f in sorted((CASE/"source").glob("*.h5ad.var.tsv"))]
    assert len(var) == 2 and var[0] == var[1] and len(set(var[0])) == len(var[0]) == 25718
    left = [full_sc_by_key[k] for k in var[0]]
    right = list(Axis([r for r in spatial if not r[1].startswith("MT-")]).by_key.values())
    # Per-feature count vectors retain the order of each deposited full spatial axis.
    n = len(spatial)
    totals = np.zeros(n, dtype=np.int64)
    slides = {}
    for row in manifest:
        if row["accession_tail"] == "530":
            with np.load(CASE/"source"/(row["file"]+".npz"), allow_pickle=False) as data:
                count = data["total"]
                assert count.shape == (n,) and np.all(count >= 0) and np.issubdtype(count.dtype,np.integer)
                assert int(count.sum(dtype=np.int64)) == int(row["total_UMI"])
                totals += count.astype(np.int64)
                slides[row["file"]] = count
    amount = {r[0]: int(v) for r, v in zip(spatial, totals)}
    assert len(amount) == n
    source_duplicates = repeated(full_sc) | repeated(full_spatial)
    metrics = analyze(left, right, source_duplicates, amount)
    write_table(out, "lung_arms.tsv", metrics)
    # Independent validation against previous AnnData-based reconstruction and its per-key count annotation.
    left_keys = {r.key: r for r in left}
    right_keys = {r.key: r for r in right}
    local_bad = {k: (left_keys[k].id, right_keys[k].id, amount[right_keys[k].id])
                 for k in left_keys.keys() & right_keys.keys() if left_keys[k].id != right_keys[k].id}
    original = {r["input_key"]: (r["sc_source_id"], r["vis_id"], sum(int(r[x+"_UMI"]) for x in ("Tumour","Background","Healthy")))
                for r in records(CASE/"discordant_model_keys.tsv")}
    assert local_bad == original and len(original)==221
    left_ids, right_ids = ({r.id: r for r in side} for side in (left,right))
    shared = left_ids.keys() & right_ids.keys()
    disagree = {i for i in shared if left_ids[i].key != right_ids[i].key}
    dropped = {i for i in shared if left_ids[i].key == right_ids[i].key and
               (left_ids[i].symbol in source_duplicates or right_ids[i].symbol in source_duplicates)}
    assert (len(disagree), sum(amount[i] for i in disagree)) == (metrics[3]["correct_ids_gained_vs_local"],metrics[3]["spatial_umi_correct_gained_vs_local"])
    assert (len(dropped), sum(amount[i] for i in dropped)) == (metrics[1]["correct_ids_lost_vs_local"],metrics[1]["spatial_umi_correct_lost_vs_local"])
    delta = [dict(comparison=event, stable_id=i, reference_symbol=left_ids[i].symbol, spatial_symbol=right_ids[i].symbol,
                  reference_key=left_ids[i].key, spatial_key=right_ids[i].key, spatial_umi=amount[i])
             for event, ids in (("local_correct_removed_by_source_drop",dropped),("shared_record_missed_by_local",disagree)) for i in sorted(ids)]
    write_table(out, "lung_changed_correct_pairs.tsv", delta)
    # A wrong apparent match need not have *either* literal ID available in the other axis.
    wrong_groups = defaultdict(list)
    for key, (ref_id, spatial_id, _) in local_bad.items():
        wrong_groups[(ref_id in right_ids,spatial_id in left_ids)].append(key)
    write_table(out,"lung_wrong_pair_counterparts.tsv", [dict(reference_id_present_in_spatial=int(a),
        spatial_id_present_in_reference=int(b), original_different_id_pairs=len(keys),
        spatial_umi=sum(local_bad[k][2] for k in keys),
        source_duplicate_pairs=sum(left_keys[k].symbol in source_duplicates or right_keys[k].symbol in source_duplicates for k in keys))
        for (a,b),keys in sorted(wrong_groups.items())])
    index = {r[0]: i for i,r in enumerate(spatial)}
    duplicate_bad = {s for k,(_,s,_) in local_bad.items() if left_keys[k].symbol in source_duplicates or right_keys[k].symbol in source_duplicates}
    all_bad = {s for _,s,_ in local_bad.values()}
    sets = dict(non_mt={r.id for r in right}, source_drop_removed_wrong=duplicate_bad,
                source_drop_retained_wrong=all_bad-duplicate_bad, source_drop_lost_correct=dropped,
                exact_id_gained_correct=disagree)
    indices = {name: [index[i] for i in ids] for name,ids in sets.items()}
    per_slide = [dict(slide=f, **{name+"_umi":int(c[ix].sum(dtype=np.int64)) for name,ix in indices.items()})
                 for f,c in sorted(slides.items())]
    for name,ids in sets.items():
        assert sum(r[name+"_umi"] for r in per_slide)==sum(amount[i] for i in ids)
    write_table(out,"lung_slide_impact.tsv",per_slide)
    # At the family level both sides contribute all retained member rows to a single unit.
    lf, rf = (defaultdict(set) for _ in range(2))
    for rows, fam in ((left,lf),(right,rf)):
        for r in rows:
            fam[r.symbol].add(r.id)
    groups = Counter()
    spatial_rows = Counter()
    ref_rows = Counter()
    umi = Counter()
    for s in lf.keys() & rf.keys():
        cardinality = "singleton" if len(lf[s]) == len(rf[s]) == 1 else "multirow_in_either"
        group = ("equal" if lf[s] == rf[s] else "different")+"_"+cardinality
        groups[group] += 1
        ref_rows[group] += len(lf[s])
        spatial_rows[group] += len(rf[s])
        umi[group] += sum(amount[i] for i in rf[s])
    write_table(out, "lung_pooled_symbols.tsv", [dict(record_membership=k, shared_symbols=groups[k], reference_rows=ref_rows[k],
               spatial_rows=spatial_rows[k], spatial_umi=umi[k], all_shared_symbols=sum(groups.values())) for k in sorted(groups)])
    return dict(raw_spatial_umi=int(totals.sum()), spatial_libraries=40, reference_libraries=81,
                raw_spatial_rows=n, raw_reference_rows=len(sc), ann_data_mismatch_rows_and_umi_verified=len(original),
                raw_reference_id_syntax="literal IDs as archived; no version-suffix normalization")


def panel(out):
    sources = records(PANEL/"panel_sources.tsv")
    freq = Counter(r["ordered_id_symbol_sha256"] for r in sources)
    axes = {}
    for digest in sorted(freq):
        path = PANEL/"panel_axes"/(digest+".tsv.gz")
        with gzip.open(path,"rb") as f:
            assert hashlib.sha256(f.read()).hexdigest() == digest
        axes[digest] = list(load(path,False,"-").by_key.values())
    summary = []
    for a, b in itertools.combinations_with_replacement(axes,2):
        left, right = axes[a], axes[b]
        master_keys = None
        if a != b:
            subset, master = sorted((left,right),key=len)
            subset_pairs = [(r.id,r.symbol) for r in subset]
            subset_ids = {r.id for r in subset}
            ordered_restriction = [(r.id,r.symbol) for r in master if r.id in subset_ids]
            if ordered_restriction == subset_pairs:
                master_keys = {r.id:r.key for r in master}
        count = freq[a]*freq[b] if a != b else freq[a]*(freq[a]-1)//2
        for row in analyze(left,right,repeated(left)|repeated(right),inherited=master_keys):
            if row["arm"] != "drop_final_duplicates":  # full, unfiltered axes here
                summary.append(dict(left_axis=a[:8],right_axis=b[:8],manifest_pairs=count,**row))
    write_table(out, "panel_arms.tsv", summary)
    local_results = {frozenset((r["left_axis"],r["right_axis"])):r for r in summary if r["arm"]=="local_keys"}
    archived = records(PANEL/"panel_pairs.tsv")
    for r in archived:
        ours = local_results[frozenset((r["left_axis"][:8],r["right_axis"][:8]))]
        for before, after in (("shared_keys","matched_pairs"),("shared_ids","exact_shared_ids"),
                              ("different_ids_duplicate","wrong_source_duplicate"),("different_ids_singleton","wrong_source_singleton"),
                              ("shared_id_different_keys","shared_ids_missed")):
            assert int(r[before]) == ours[after]
    return dict(tables=len(sources), axes=len(axes), different_axis_comparisons=10, self_axis_controls=5,
                manifest_pairs=sum(r["manifest_pairs"] for r in summary if r["arm"]=="local_keys"),
                checked_archived_panel_comparisons=len(archived))


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path,default=HERE)
    out=parser.parse_args().output_dir
    out.mkdir(exist_ok=True,parents=True)
    checks={"lung":lung(out),"panel":panel(out)}
    (out/"validation.json").write_text(json.dumps(checks,indent=2)+"\n",encoding="utf-8")
