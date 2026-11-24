"""Audit two gene axes without modifying data. Standard library only, Python 3.10+."""

import argparse
import csv
import gzip
import json
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Row:
    id: str
    symbol: str
    key: str


def assign_keys(symbols, separator):
    # AnnData retains original strings (including existing suffix-like symbols).
    reserved = set(symbols)
    seen = set()
    last = Counter()
    for name in symbols:
        if name not in seen:
            seen.add(name)
            yield name
            continue
        while True:
            last[name] += 1
            candidate = f"{name}{separator}{last[name]}"
            if candidate not in reserved:
                reserved.add(candidate)
                yield candidate
                break


class Axis:
    def __init__(self, entries, separator="-"):
        if not entries:
            raise ValueError("No gene records supplied")
        keys = [r[2] for r in entries]
        if all(k is None for k in keys):
            keys = assign_keys([r[1] for r in entries], separator)
        self.by_key = {}
        self.by_id = defaultdict(list)
        self.by_symbol = defaultdict(set)
        self.multiplicity = Counter()
        for (fid, symbol, _), key in zip(entries, keys):
            if not all((fid, symbol, key)) or key in self.by_key:
                raise ValueError(f"Missing ID, symbol or key, or duplicate final key: {(fid, symbol, key)}")
            row = Row(fid, symbol, key)
            self.by_key[key] = row
            self.by_id[fid].append(row)
            self.by_symbol[symbol].add(fid)
            self.multiplicity[symbol] += 1

    def use_source_groups(self, source):
        """Retain analyzed rows/keys but classify duplicates on their complete source axis."""
        for row in self.by_key.values():
            if row.id not in source.by_symbol.get(row.symbol, ()):
                raise ValueError(f"Analyzed feature absent from supplied provenance: {row.id}, {row.symbol}")
        self.by_symbol = source.by_symbol
        self.multiplicity = source.multiplicity


def parse(handle, keyed=False):
    rows = []
    reader = csv.reader(handle, delimiter="\t")
    if keyed and next(reader, None) != ["stable_id", "symbol", "final_key"]:
        raise ValueError("Expected header stable_id, symbol, final_key (tab separated)")
    for number, fields in enumerate(reader, 2 if keyed else 1):
        if keyed and len(fields) == 3:
            rows.append(tuple(fields))
        elif not keyed and len(fields) in (2, 3):
            if len(fields) == 3 and fields[2] != "Gene Expression":
                raise ValueError("Raw mode expects a gene-only axis. For mixed modalities, export the gene rows "
                                 "with their actual final keys and use keyed mode.")
            rows.append((fields[0], fields[1], None))
        else:
            raise ValueError(f"Unexpected TSV columns on line {number}: {len(fields)}")
    return rows


def load(path, keyed, separator):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8-sig", newline="") as handle:
        return Axis(parse(handle, keyed), separator)


EVENT_FIELDS = ("event", "class", "left_id", "right_id", "left_symbol", "right_symbol",
                "left_key", "right_key", "left_multiplicity", "right_multiplicity",
                "cross_table_family_id_overlap", "left_key_absent_right", "right_key_absent_left")


def analyze(left, right):
    shared_keys = left.by_key.keys() & right.by_key.keys()
    shared_ids = left.by_id.keys() & right.by_id.keys()
    summary = dict(shared_keys=len(shared_keys), shared_ids=len(shared_ids),
                   different_ids_duplicate=0, different_ids_singleton=0, duplicate_with_family_id_overlap=0,
                   shared_ids_nonunique_in_either_axis=0, shared_id_different_keys=0,
                   shared_id_changed_symbol_key_preserved=0,
                   singly_occurring_renamed_ids=0, singleton_renames_both_keys_absent_opposite=0)
    events = []

    def add(event, kind, a, b, overlap):
        events.append((event, kind, a.id, b.id, a.symbol, b.symbol, a.key, b.key,
                       left.multiplicity[a.symbol], right.multiplicity[b.symbol], overlap,
                       a.key not in right.by_key, b.key not in left.by_key))

    for key in sorted(shared_keys):
        a, b = left.by_key[key], right.by_key[key]
        if a.id == b.id:
            continue
        duplicate = left.multiplicity[a.symbol] > 1 or right.multiplicity[b.symbol] > 1
        overlap = a.id in right.by_symbol.get(a.symbol, ()) or b.id in left.by_symbol.get(b.symbol, ())
        summary["different_ids_duplicate" if duplicate else "different_ids_singleton"] += 1
        summary["duplicate_with_family_id_overlap"] += duplicate and overlap
        add("same_key_different_id", "duplicate_associated" if duplicate else "singleton_on_both", a, b, overlap)
    for fid in sorted(shared_ids):
        if len(left.by_id[fid]) != 1 or len(right.by_id[fid]) != 1:
            summary["shared_ids_nonunique_in_either_axis"] += 1
            continue
        a, b = left.by_id[fid][0], right.by_id[fid][0]
        if a.key == b.key:
            if a.symbol != b.symbol:
                summary["shared_id_changed_symbol_key_preserved"] += 1
                add("same_id_changed_symbol_key_preserved", "annotation_change_only", a, b, "")
            continue
        summary["shared_id_different_keys"] += 1
        unique_rename = (left.multiplicity[a.symbol] == right.multiplicity[b.symbol] == 1
                         and a.symbol != b.symbol)
        summary["singly_occurring_renamed_ids"] += unique_rename
        summary["singleton_renames_both_keys_absent_opposite"] += (
            unique_rename and a.key not in right.by_key and b.key not in left.by_key)
        add("same_id_different_keys", "singleton_rename" if unique_rename else "other_key_change", a, b, "")
    return summary, events


def main():
    arg = argparse.ArgumentParser(description=__doc__)
    arg.add_argument("left", type=Path)
    arg.add_argument("right", type=Path)
    arg.add_argument("--left-keyed", action="store_true", help="Input supplies actual final keys")
    arg.add_argument("--right-keyed", action="store_true", help="Input supplies actual final keys")
    arg.add_argument("--left-source-features", type=Path, help="Full original gene table for left duplicate provenance")
    arg.add_argument("--right-source-features", type=Path, help="Full original gene table for right duplicate provenance")
    arg.add_argument("--separator", default="-", help="Separator when generating AnnData-style keys")
    arg.add_argument("--events", type=Path, required=True, help="Output discrepancy TSV")
    args = arg.parse_args()
    left = load(args.left, args.left_keyed, args.separator)
    right = load(args.right, args.right_keyed, args.separator)
    for axis, source in ((left, args.left_source_features), (right, args.right_source_features)):
        if source:
            axis.use_source_groups(load(source, False, args.separator))
    summary, events = analyze(left, right)
    inputs = (args.left, args.right, args.left_source_features, args.right_source_features)
    if args.events.resolve() in {p.resolve() for p in inputs if p}:
        arg.error("Use a separate path for the events file")
    with args.events.open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(EVENT_FIELDS)
        writer.writerows(events)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
