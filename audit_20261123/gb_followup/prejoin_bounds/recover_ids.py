"""Bound a row's source identity from names, a known source axis and optional order."""

import argparse
from bisect import bisect_left, bisect_right
from collections import defaultdict
import csv
import gzip
from pathlib import Path
import re


def read_text(path):
    opener = gzip.open if path.suffix == ".gz" else open
    return opener(path, "rt", encoding="utf-8-sig", newline="")


class Reference:
    def __init__(self, records):
        self.ids = []
        self.positions = defaultdict(list)
        for row, (identifier, symbol) in enumerate(records):
            if not identifier or not symbol:
                raise ValueError("Empty source identifier or symbol")
            self.ids.append(identifier)
            self.positions[symbol].append(row)
        if not self.ids or len(set(self.ids)) != len(self.ids):
            raise ValueError("Source must contain nonempty, one-to-one identifiers")

    def possible_rows(self, label):
        # Preserve real source names ending in e.g. "-1". A generated suffix
        # requires two instances of an identical original name on this source.
        possible = self.positions.get(label, [])
        suffix = re.fullmatch(r"(.+)-[1-9][0-9]*", label)
        if suffix:
            base = self.positions.get(suffix[1], [])
            if len(base) > 1:
                return sorted([*possible, *base])
        return possible

    def candidates(self, labels, preserve_order=False):
        """Outer sets of compatible IDs. Inputs need no counts or retained IDs."""
        if not labels or len(set(labels)) != len(labels):
            raise ValueError("Supply the complete ordered list of nonempty, unique final keys")
        rows = [self.possible_rows(key) for key in labels]
        for key, choices in zip(labels, rows):
            if not choices:
                raise ValueError(f"No compatible source symbol for {key!r}")
        if preserve_order:
            # Earliest possible end of each prefix. Each retained source row
            # must come after this prefix and before a feasible suffix.
            earliest = []
            previous = -1
            for key, choices in zip(labels, rows):
                index = bisect_right(choices, previous)
                if index == len(choices):
                    raise ValueError(f"No order-preserving source embedding at {key!r}")
                previous = choices[index]
                earliest.append(previous)
            following = len(self.ids)
            for j in range(len(rows) - 1, -1, -1):
                choices = rows[j]
                end = bisect_left(choices, following)
                start = bisect_right(choices, earliest[j-1] if j else -1)
                rows[j] = choices[start:end]
                following = choices[end-1]
        return {label: tuple(self.ids[i] for i in choices) for label, choices in zip(labels, rows)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path, help="Headerless gene-only ID, symbol [, Gene Expression] TSV")
    parser.add_argument("keys", type=Path, help="One final key per line, no header")
    parser.add_argument("output", type=Path, help="Output TSV, one row per candidate ID")
    parser.add_argument("--order-preserved", action="store_true", help="Historical source order is known to be preserved")
    args = parser.parse_args()
    if args.output.resolve() in {args.source.resolve(), args.keys.resolve()}:
        parser.error("Choose an output path distinct from both inputs")
    with read_text(args.source) as f:
        records = []
        for row in csv.reader(f, delimiter="\t"):
            if len(row) != 2 and not (len(row) == 3 and row[2] == "Gene Expression"):
                raise ValueError("Supply a headerless, single-modality gene table")
            records.append(row[:2])
    with read_text(args.keys) as f:
        keys = [line.rstrip("\r\n") for line in f]
    candidates = Reference(records).candidates(keys, args.order_preserved)
    with args.output.open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(["final_key", "candidate_count", "source_id"])
        for label, ids in candidates.items():
            writer.writerows((label, len(ids), identifier) for identifier in ids)


if __name__ == "__main__":
    main()
