"""Replay masked-ID tests on all selected axes, plus exhaustive small controls."""

import argparse
import csv
import gzip
from itertools import combinations, product
import json
from pathlib import Path
import sys
import tempfile

from recover_ids import Reference, read_text

HERE = Path(__file__).resolve().parent
F = HERE.parent
sys.path.insert(0, str(F / "systematic_mechanism"))
from audit_join import assign_keys


def table(path):
    with read_text(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def raw(path):
    with read_text(path) as f:
        return [tuple(row[:2]) for row in csv.reader(f, delimiter="\t")]


def generate_axis(records):
    ids = [identifier for identifier, _ in records]
    keys = list(assign_keys([sym for _, sym in records], "-"))
    return keys, dict(zip(keys, ids))


def control_checks():
    # Exhaustively compare candidate bounds against all actual ordered subsets
    # for both naming orders, including literal/generated suffix collisions.
    checked = 0
    for syms in product(("X", "X-1", "Y"), repeat=4):
        records = [(str(i), x) for i, x in enumerate(syms)]
        reference = Reference(records)
        inherited, _ = generate_axis(records)
        for size in range(1, 5):
            for keep in combinations(range(4), size):
                subset = [records[i] for i in keep]
                filtered, _ = generate_axis(subset)
                for keys in (filtered, [inherited[i] for i in keep]):
                    bounds = reference.candidates(keys, preserve_order=True)
                    assert all(str(i) in bounds[key] for i, key in zip(keep, keys))
                    checked += 1
    adjacent = Reference([("1", "X"), ("2", "X")]).candidates(["X"], True)
    assert adjacent["X"] == ("1", "2")
    # Deliberate reversal with two unique anchors.
    try:
        Reference([("1", "A"), ("2", "B")]).candidates(["B", "A"], True)
    except ValueError as error:
        assert "No order-preserving" in str(error)
    else:
        raise AssertionError("Reversed unique anchors accepted")
    return checked


def save_tsv(path, records):
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(records[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)


def run(out):
    tests = control_checks()
    axis_rows, ambigs, pairs = [], [], []

    def evaluate(reference, observed, dataset, sample, threshold):
        keys, gold = generate_axis(observed)
        # The reconstruction never receives the gold dictionary or the retained rows.
        un = reference.candidates(keys)
        ordered = reference.candidates(keys, True)
        n_ambig = n_left = verified = 0
        for key in keys:
            choices = ordered[key]
            assert set(choices) <= set(un[key]) and gold[key] in choices
            if len(un[key]) > 1:
                n_ambig += 1
            if len(choices) > 1:
                ambigs.append(dict(dataset=dataset, sample=sample, threshold=threshold, key=key,
                                   candidate_ids="|".join(choices), gold_id=gold[key]))
                n_left += 1
            else:
                verified += int(choices[0] == gold[key])
        axis_rows.append(dict(dataset=dataset, sample=sample, threshold=threshold,
                              rows=len(keys), ambiguous_without_order=n_ambig, ambiguous_with_order=n_left,
                              assigned_correctly=verified, wrong_singletons=len(keys)-n_left-verified))
        return gold, ordered, un

    src = F / "systematic_mechanism" / "panel_axes"
    reference = Reference(raw(next(src.glob("dd4b9*.tsv.gz"))))
    observe = raw(next(src.glob("1400d*.tsv.gz")))
    ff_gold, ff_order, ff_unordered = evaluate(reference, observe, "probe", "FFPE", "panel")
    ff_key = "ARMCX5-GPRASP2"
    ff = dict(key=ff_key, unordered_candidates=ff_unordered[ff_key], ordered_candidates=ff_order[ff_key],
              true_id=ff_gold[ff_key])

    cortex = F / "external_order_replication"
    cortex_rows = table(cortex / "gene_totals.tsv.gz")
    cortex_raw = [(r["id"], r["symbol"]) for r in cortex_rows]
    cortex_ref = Reference(cortex_raw)
    for info in table(cortex / "sections.tsv"):
        key = info["section"]
        evaluate(cortex_ref, [item for item,r in zip(cortex_raw,cortex_rows) if int(r[key])>0],
                 "cortex", key, 1)

    blood = F / "prospective_holdout"
    blood_rows = table(blood / "axis_counts.tsv.gz")
    blood_raw = [(r["id"], r["symbol"]) for r in blood_rows]
    blood_ref = Reference(blood_raw)
    samples = table(blood / "samples.tsv")
    for threshold in (1, 3, 10):
        within = {}
        for sample in samples:
            label = sample["sample"]
            survivors = [item for item, row in zip(blood_raw, blood_rows)
                         if int(row[label+"_cells"]) >= threshold]
            gold, ordered, _ = evaluate(blood_ref, survivors, "PBMC", label, threshold)
            if threshold == 1:
                within[label] = (gold, ordered)
        if threshold != 1:
            continue
        for condition in ("ctrl", "stim"):
            labels = [r["sample"] for r in samples if r["condition"] == condition]
            for left, right in combinations(labels, 2):
                a, ca = within[left]
                b, cb = within[right]
                same = different = unresolved = errors = 0
                misses = []
                common = a.keys() & b.keys()
                for key in common:
                    la, rb = ca[key], cb[key]
                    if la == rb and len(la) == 1:
                        same += 1
                        errors += int(a[key] != b[key])
                    elif set(la).isdisjoint(rb):
                        different += 1
                        errors += int(a[key] == b[key])
                    else:
                        unresolved += 1
                        if a[key] != b[key]:
                            misses.append(key)
                assert not errors and same + different + unresolved == len(common)
                pairs.append(dict(condition=condition, left=left, right=right, shared=len(common),
                                  certified_same=same, certified_different=different, unresolved=unresolved,
                                  wrong_certifications=errors, actual_differences=different+len(misses),
                                  unresolved_actual_keys="|".join(sorted(misses))))

    # Real sorted final axes fail the ordering precondition. No IDs inferred
    # from row order are used for these publications.
    lung = F / "fifth_case"
    sources = {r["axis_id"] for r in table(lung/"input_axes.tsv") if r["accession_tail"] == "526"}
    assert len(sources) == 1
    lung_ref = Reference(raw(lung / "source" / (sources.pop()+".tsv.gz")))
    lung_out = {}
    for axis in sorted((lung/"source").glob("*.h5ad.var.tsv")):
        keys = [r["_index"] for r in table(axis)]
        try:
            lung_ref.candidates(keys, True)
        except ValueError as exc:
            assert "No order-preserving" in str(exc)
            lung_out[axis.name] = str(exc)
        else:
            raise AssertionError("Unexpected compatible source order in lung archive")
    assert len(lung_out) == 2
    save_tsv(out/"axes.tsv", axis_rows)
    save_tsv(out/"pairs.tsv", pairs)
    save_tsv(out/"ambiguous.tsv", ambigs)
    result = dict(masked_axes=len(axis_rows), small_axis_tests=tests, probe_sentinel=ff,
                  incompatible_lung_axes=lung_out, axis_aggregates={}, pbmc_threshold1={})
    for dataset in ("probe","cortex","PBMC"):
        subset = [r for r in axis_rows if r["dataset"] == dataset]
        result["axis_aggregates"][dataset] = {k:sum(r[k] for r in subset) for k in
            ("rows","ambiguous_without_order","ambiguous_with_order","assigned_correctly","wrong_singletons")}
    for condition in ("ctrl","stim"):
        subset = [r for r in pairs if r["condition"] == condition]
        result["pbmc_threshold1"][condition] = {k:sum(r[k] for r in subset) for k in
            ("shared","certified_different","unresolved","actual_differences","wrong_certifications")}
    (out/"summary.json").write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true", help="Regenerate the frozen result files")
    args = parser.parse_args()
    if args.write:
        run(HERE)
    else:
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)
            run(path)
            for file in path.iterdir():
                assert file.read_text(encoding="utf8") == (HERE/file.name).read_text(encoding="utf8"), file.name
        print("Masked-axis counts, candidates, pair decisions and failure controls reproduced.")
