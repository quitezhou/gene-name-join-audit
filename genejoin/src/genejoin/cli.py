"""Command-line interface.

Exit status for ``check``: 0 no discordance, 1 discordances found, 2 input error.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
from contextlib import contextmanager
from pathlib import Path
from tempfile import NamedTemporaryFile

from . import __version__
from .core import EVENT_FIELDS, Reference, audit
from .io import load_axis, read_feature_table, read_keys, resolve_input


@contextmanager
def _output_table(path, inputs):
    output = Path(path)
    for source in (resolve_input(p) for p in inputs if p):
        if output.resolve() == source.resolve() or (output.exists() and output.samefile(source)):
            raise ValueError("Choose an output path distinct from every input file")
    temp = NamedTemporaryFile(mode="w", encoding="utf-8", newline="", dir=output.parent, delete=False)
    try:
        with temp:
            yield temp
        os.replace(temp.name, output)
    finally:
        Path(temp.name).unlink(missing_ok=True)


def _check(args) -> int:
    opts = dict(separator=args.separator, mixed_modalities=args.mixed_modalities,
                id_column=args.id_column, symbol_column=args.symbol_column)
    left = load_axis(args.left, args.left_format, **opts)
    right = load_axis(args.right, args.right_format, **opts)
    for axis, source in ((left, args.left_source), (right, args.right_source)):
        if source:
            axis.use_source_groups(load_axis(source, "features", **opts))
    report = audit(left, right)
    if args.events:
        with _output_table(args.events, (args.left, args.right, args.left_source, args.right_source)) as f:
            w = csv.DictWriter(f, fieldnames=EVENT_FIELDS, delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(report.events)
    if args.json:
        print(json.dumps(report.to_dict(), indent=2))
    else:
        s = report.summary
        print(f"shared names: {s['shared_keys']:,}   shared IDs: {s['shared_ids']:,}")
        print(f"same name, different ID: {report.same_name_different_id:,} "
              f"({s['different_ids_duplicate']:,} duplicate-associated, {s['different_ids_singleton']:,} singleton, "
              f"{s.get('different_ids_unclassified', 0):,} unclassified)")
        print(f"same ID, different name: {report.same_id_different_name:,} "
              f"({s['singly_occurring_renamed_ids']:,} singleton symbol changes)")
        if s["shared_ids_nonunique_in_either_axis"]:
            print(f"shared IDs occurring more than once on a side (not checked): "
                  f"{s['shared_ids_nonunique_in_either_axis']:,}")
        print("OK: no discordance detected among comparable rows." if report.ok else
              "DISCORDANT: review the mismatched records before joining.")
    return 0 if report.ok else 1


def _recover(args) -> int:
    source = read_feature_table(args.source, args.mixed_modalities)
    keys = read_keys(args.keys)
    result = Reference(source).candidates(keys, args.order_preserved)
    with _output_table(args.output, (args.source, args.keys)) as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["final_key", "candidate_count", "source_id"])
        for label, ids in result.items():
            w.writerows((label, len(ids), i) for i in ids)
    resolved = sum(len(v) == 1 for v in result.values())
    print(f"{len(result):,} keys: {resolved:,} resolve to one source ID, {len(result) - resolved:,} remain ambiguous")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="genejoin", description="Check gene-name joins between datasets.")
    p.add_argument("--version", action="version", version=f"genejoin {__version__}")
    sub = p.add_subparsers(dest="command", required=True)

    c = sub.add_parser("check", help="Audit a name join between two datasets",
                       description="Compare each shared name with its stable IDs, and each shared ID with its "
                                   "names. Inputs: 10x feature table or folder (keys generated in file order, "
                                   "AnnData style), keyed TSV (stable_id, symbol, final_key) or .h5ad.")
    c.add_argument("left")
    c.add_argument("right")
    for side in ("left", "right"):
        c.add_argument(f"--{side}-format", choices=["auto", "features", "keyed", "h5ad"], default="auto")
        c.add_argument(f"--{side}-source", help="Complete unfiltered feature table this side came from "
                                                "(improves duplicate-family classification)")
    c.add_argument("--events", help="Write one row per discordance to this TSV")
    c.add_argument("--json", action="store_true", help="Print the summary as JSON")
    c.add_argument("--separator", default="-", help="Suffix separator when generating keys (default '-')")
    c.add_argument("--mixed-modalities", choices=["refuse", "genes-only"], default="refuse")
    c.add_argument("--id-column", help=".h5ad: column of .var holding stable IDs")
    c.add_argument("--symbol-column", help=".h5ad: column of .var holding original symbols")
    c.set_defaults(func=_check)

    r = sub.add_parser("recover", help="Bound source IDs for keys of an object that lost its ID column",
                       description="Given the original feature table and the final ordered keys, list the source "
                                   "IDs each key could denote. Use --order-preserved only when provenance shows "
                                   "that rows kept their source order.")
    r.add_argument("source", help="Original headerless feature table (ID, symbol [, type])")
    r.add_argument("keys", help="Final keys, one per line, in their actual order")
    r.add_argument("output", help="Output TSV: final_key, candidate_count, source_id")
    r.add_argument("--order-preserved", action="store_true")
    r.add_argument("--mixed-modalities", choices=["refuse", "genes-only"], default="refuse")
    r.set_defaults(func=_recover)
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return args.func(args)
    except (ValueError, OSError, EOFError, csv.Error) as exc:
        print(f"genejoin: error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
