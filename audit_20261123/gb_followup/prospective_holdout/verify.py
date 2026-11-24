"""Replay locked rank predictions against archived AnnData outputs (stdlib only)."""
import argparse
import csv
import gzip
import hashlib
import io
import itertools
import json
import re
import statistics
import warnings
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
THRESHOLDS = (0, 1, 3, 10)


def tsv_bytes(header, rows):
    stream = io.StringIO(newline="")
    writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
    writer.writerow(header)
    writer.writerows(rows)
    return stream.getvalue().encode("utf-8")


def table(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as fin:
        return list(csv.DictReader(fin, delimiter="\t"))


def rank_names(axis):
    """Predict the name of the kth retained member after reserving native labels."""
    families = defaultdict(list)
    for identifier, name in axis:
        families[name].append(identifier)
    result = {}
    for base, members in families.items():
        result[members[0]] = base
        candidate = 1
        for identifier in members[1:]:
            while f"{base}-{candidate}" in families:
                candidate += 1
            result[identifier] = f"{base}-{candidate}"
            candidate += 1
    return result


def predictable(left, right):
    right_at = {key: identifier for identifier, key in right.items()}
    collisions = set()
    splits = set()
    key_count, id_count = 0, 0
    for identifier, key in left.items():
        if key in right_at:
            key_count += 1
            if right_at[key] != identifier:
                collisions.add(("mismatch", key, identifier, right_at[key]))
        if identifier in right:
            id_count += 1
            if right[identifier] != key:
                splits.add(("split", identifier, key, right[identifier]))
    return collisions | splits, key_count, id_count


def name_hash(mapping):
    rows = sorted(mapping.items())
    return hashlib.sha256(tsv_bytes(("id", "key"), rows)).hexdigest()


def generate_official(retained, pairs):
    import anndata
    import numpy as np
    import pandas as pd
    if anndata.__version__ != "0.11.4":
        raise ValueError(f"This replay specifies AnnData 0.11.4; got {anndata.__version__}")
    official = {}
    hashes = {}
    events = []
    for (threshold, sample), axis in retained.items():
        ids, names = zip(*axis)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            obj = anndata.AnnData(np.empty((0, len(ids)), dtype=np.uint8),
                                  var=pd.DataFrame({"stable": ids}, index=names))
            obj.var_names_make_unique()
        frame = obj.var.set_axis(obj.var_names, axis=0).reset_index()
        frame.columns = ["key", "stable"]
        assert frame.key.is_unique and frame.stable.is_unique
        official[threshold, sample] = frame
        hashes[f"{threshold}/{sample}"] = name_hash(dict(zip(frame.stable, frame.key)))
    for threshold in THRESHOLDS:
        for group, left, right in pairs:
            a, b = official[threshold, left], official[threshold, right]
            prefix = (threshold, group, left, right)
            joined_names = a.merge(b, on="key", validate="1:1")
            for key, id_a, id_b in joined_names.itertuples(index=False, name=None):
                if id_a != id_b:
                    events.append(prefix + ("mismatch", key, id_a, id_b))
            joined_ids = a.merge(b, on="stable", validate="1:1")
            assert len(joined_ids) == len(set(a.stable) & set(b.stable))
            for key_a, identifier, key_b in joined_ids.itertuples(index=False, name=None):
                if key_a != key_b:
                    events.append(prefix + ("split", identifier, key_a, key_b))
    return sorted(events), hashes


def main(use_official, record):
    if record and not use_official:
        raise ValueError("Recording requires regeneration using AnnData")
    axis = table(HERE / "axis_counts.tsv.gz")
    for field in ("id", "symbol"):
        assert all(row[field] for row in axis)
    ids = [row["id"] for row in axis]
    assert len(ids) == len(set(ids))
    original = [row["symbol"] for row in axis]
    exact = Counter(original)
    folded = defaultdict(set)
    for name in exact:
        folded[name.casefold()].add(name)
    samples = table(HERE / "samples.tsv")
    names = [row["sample"] for row in samples]
    donors = {r["sample"]: r["donor"] for r in samples}
    by_group = {condition: [r["sample"] for r in samples if r["condition"] == condition]
                for condition in ("ctrl", "stim")}
    pairs = [(group, a, b) for group, members in by_group.items()
             for a, b in itertools.combinations(members, 2)]
    pairs += [("paired", a, b) for a in by_group["ctrl"] for b in by_group["stim"]
              if donors[a] == donors[b]]
    umi = {sample: {r["id"]: int(r[sample + "_umis"]) for r in axis} for sample in names}
    detects = {sample: [int(r[sample + "_cells"]) for r in axis] for sample in names}
    retained = {}
    predicted = {}
    hashes = {}
    events, pair_rows, burdens, totals = [], [], [], []
    for threshold in THRESHOLDS:
        all_affected = defaultdict(set)
        for sample in names:
            kept = [(i, s) for i, s, n in zip(ids, original, detects[sample]) if n >= threshold]
            retained[threshold, sample] = kept
            named = rank_names(kept)
            assert len(set(named.values())) == len(kept)
            predicted[sample] = named
            hashes[f"{threshold}/{sample}"] = name_hash(named)
        for group, a, b in pairs:
            result, key_count, id_count = predictable(predicted[a], predicted[b])
            prefix = (threshold, group, a, b)
            events.extend(prefix + item for item in result)
            error_counts = Counter(event[0] for event in result)
            pair_rows.append(prefix + (key_count, id_count, error_counts["mismatch"], error_counts["split"]))
            if threshold == 0:
                assert not result
            for kind, index, x, y in result:
                left_id, right_id = (x, y) if kind == "mismatch" else (index, index)
                all_affected[group, a, kind].add(left_id)
                all_affected[group, b, kind].add(right_id)
        for group in ("ctrl", "stim", "paired"):
            members = by_group[group] if group in by_group else names
            for sample in members:
                mismatch = all_affected[group, sample, "mismatch"]
                splits = all_affected[group, sample, "split"]
                for metric, affected in (("mismatch", mismatch), ("split", splits), ("union", mismatch | splits)):
                    retained_umi = sum(umi[sample][identifier] for identifier in predicted[sample])
                    affected_umi = sum(umi[sample][identifier] for identifier in affected)
                    burdens.append((threshold, group, sample, donors[sample], metric, len(affected),
                                    len(predicted[sample]), affected_umi, retained_umi))
            lines = [r for r in pair_rows if r[:2] == (threshold, group)]
            group_events = [r for r in events if r[:2] == (threshold, group)]
            affected_donors = {donors[s] for s in members
                               if all_affected[group, s, "mismatch"] or all_affected[group, s, "split"]}
            selected = [r for r in burdens if r[:2] == (threshold, group) and r[4] == "union"]
            donor_umis = defaultdict(lambda: [0, 0])
            for row in selected:
                donor_umis[row[3]][0] += row[-2]
                donor_umis[row[3]][1] += row[-1]
            fractions_ppm = [1e6 * n / d for n, d in donor_umis.values()]
            rates = [round(v, 4) for v in (min(fractions_ppm), statistics.median(fractions_ppm), max(fractions_ppm))]
            totals.append((threshold, group, len(lines), sum(r[6] > 0 for r in lines),
                           sum(r[7] > 0 for r in lines), sum(r[6] + r[7] > 0 for r in lines),
                           len(affected_donors), len({donors[s] for s in members}),
                           len({r[5] for r in group_events if r[4] == "mismatch"}),
                           len({r[5] for r in group_events if r[4] == "split"}),
                           sum(r[-2] for r in selected), sum(r[-1] for r in selected), *rates))
    events.sort()
    events_header = ("threshold", "comparison", "left", "right", "event", "index", "left_value", "right_value")
    event_bytes = tsv_bytes(events_header, events)
    if use_official:
        official_events, official_hashes = generate_official(retained, pairs)
        actual = set(official_events)
        expected = set(events)
        if actual != expected:
            raise ValueError(f"Unexpected events {len(actual-expected)}, missing {len(expected-actual)}")
        if hashes != official_hashes:
            raise ValueError("Full official stable-ID to local-name map differs")
    events_file = HERE / "official_events.tsv.gz"
    hash_file = HERE / "official_name_hashes.json"
    if record and not events_file.exists():
        with events_file.open("wb") as raw:
            with gzip.GzipFile(fileobj=raw, filename="", mode="wb", mtime=0) as compressed:
                compressed.write(tsv_bytes(events_header, official_events))
        hash_file.write_text(json.dumps(official_hashes, indent=2) + "\n", encoding="utf-8")
    assert gzip.decompress(events_file.read_bytes()) == event_bytes
    assert json.loads(hash_file.read_text()) == hashes
    diagnostic = {
        "features": len(axis), "distinct_raw_symbols": len(exact),
        "raw_duplicate_symbol_families": sum(n > 1 for n in exact.values()),
        "casefold_only_multispelling_families": sum(len(s) > 1 for s in folded.values()),
        "raw_natural_suffix_spellings_with_present_base": sum(
            bool((match := re.fullmatch(r"(.+)-([1-9][0-9]*)", key)) and match[1] in exact) for key in exact),
        "donors": len(set(donors.values())), "sample_conditions": len(names),
        "singlets": sum(int(row["singlets"]) for row in samples), "pair_threshold_comparisons": len(pair_rows),
        "map_digests_compared": len(hashes), "event_records_matched": len(events),
        "source_protocol_sha256": hashlib.sha256((HERE / "protocol.md").read_bytes()).hexdigest(),
    }
    reports = {
        "pair_summary.tsv": tsv_bytes(("threshold", "comparison", "left", "right", "shared_keys",
            "shared_ids", "mismatched_keys", "split_ids"), pair_rows),
        "burden_by_sample.tsv": tsv_bytes(("threshold", "comparison", "sample", "donor", "event",
            "affected_ids", "retained_ids", "affected_UMI", "retained_UMI"), burdens),
        "summary.tsv": tsv_bytes(("threshold", "comparison", "pairs", "pairs_with_mismatch",
            "pairs_with_split", "pairs_with_either", "donors_with_either", "all_donors",
            "distinct_mismatched_keys", "distinct_split_ids", "union_affected_UMI", "retained_UMI",
            "min_donor_ppm", "median_donor_ppm", "max_donor_ppm"), totals),
        "checks.json": (json.dumps(diagnostic, indent=2) + "\n").encode(),
    }
    for filename, data in reports.items():
        path = HERE / filename
        if record and not path.exists():
            path.write_bytes(data)
        if path.read_bytes() != data:
            raise ValueError(f"Replay differs from archived {filename}")
    print(json.dumps(diagnostic))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--official", action="store_true", help="repeat reference with AnnData 0.11.4")
    parser.add_argument("--record", action="store_true", help="create missing reference files on first run")
    args = parser.parse_args()
    main(args.official, args.record)
