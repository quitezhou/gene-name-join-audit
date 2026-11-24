"""Rebuild paired gene-axis checks from five bundled axes or the original 13 tables."""
import argparse
import csv
import gzip
import hashlib
import io
import json
import sys
from collections import Counter
from itertools import combinations
from pathlib import Path

from audit_join import Axis, EVENT_FIELDS, analyze, load, parse

HERE = Path(__file__).resolve().parent


def compressed(path, text):
    with path.open("wb") as output, gzip.GzipFile(fileobj=output, filename="", mtime=0, mode="wb") as stream:
        stream.write(text.encode("utf-8"))


def write_table(path, records):
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(records[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)


def extract_manifest(path, dest):
    # Reuse the project's established reader for the three legacy tarred matrices.
    sys.path.insert(0, str(path.parent))
    from audit_feature_panel import open_text
    rows = []
    axes_dir = dest / "panel_axes"
    axes_dir.mkdir(exist_ok=True)
    fingerprints = set()
    for dataset in json.loads(path.read_text(encoding="utf-8"))["datasets"]:
        src = path.parent / dataset["local_path"]
        stream, member = open_text(src)
        with stream:
            entries = parse(stream)
        normalized = "".join(f"{i}\t{s}\n" for i,s,_ in entries)
        digest = hashlib.sha256(normalized.encode("utf-8")).hexdigest()
        if digest not in fingerprints:
            compressed(axes_dir / f"{digest}.tsv.gz", normalized)
            fingerprints.add(digest)
        rows.append(dict(label=dataset["label"], project=dataset["project"], source_url=dataset["source_url"],
                         local_path=dataset["local_path"], selected_member=member, gene_rows=len(entries),
                         ordered_id_symbol_sha256=digest))
    write_table(dest / "panel_sources.tsv", rows)


def ordered_subset_control(axes):
    nested = []
    for a, b in combinations(axes, 2):
        big, small = sorted((a, b), key=lambda k: len(axes[k].by_key), reverse=True)
        full, part = axes[big], axes[small]
        retained = [r for r in full.by_key.values() if r.id in part.by_id]
        if [(r.id, r.symbol) for r in retained] != [(r.id, r.symbol) for r in part.by_key.values()]:
            continue
        # The same retained records receive their keys at different stages.
        inherited = Axis([(r.id, r.symbol, r.key) for r in retained])
        stats, events = analyze(full, part)
        restored, fixes = analyze(full, inherited)
        assert not fixes and restored["shared_keys"] == len(retained)
        assert stats["shared_ids"] == len(retained) and stats["shared_ids_nonunique_in_either_axis"] == 0
        nested.append(dict(full_fingerprint=big, subset_fingerprint=small, full_rows=len(full.by_key),
                           retained_rows=len(retained), independent_name_matches=stats["shared_keys"],
                           independent_same_key_wrong_id=stats["different_ids_duplicate"] + stats["different_ids_singleton"],
                           inherited_name_matches=restored["shared_keys"], inherited_same_key_wrong_id=0,
                           direct_one_to_one_id_matches=stats["shared_ids"],
                           independent_discrepancies=[dict(zip(EVENT_FIELDS, e)) for e in events]))
    return nested


def rebuild(bundle, dest):
    with (bundle / "panel_sources.tsv").open(encoding="utf-8", newline="") as f:
        sources = list(csv.DictReader(f, delimiter="\t"))
    axes = {}
    for row in sources:
        digest = row["ordered_id_symbol_sha256"]
        path = bundle / "panel_axes" / f"{digest}.tsv.gz"
        if digest not in axes:
            if hashlib.sha256(gzip.decompress(path.read_bytes())).hexdigest() != digest:
                raise ValueError(f"Unexpected contents in {path}")
            axes[digest] = load(path, False, "-")
        assert len(axes[digest].by_key) == int(row["gene_rows"])

    rows = []
    results = {}
    buf = io.StringIO(newline="")
    details = csv.writer(buf, delimiter="\t", lineterminator="\n")
    details.writerow(("left_axis", "right_axis", *EVENT_FIELDS))
    for left, right in combinations(axes, 2):
        summary, events = analyze(axes[left], axes[right])
        rows.append(dict(left_axis=left, right_axis=right, **summary))
        results[frozenset((left,right))] = summary
        details.writerows((left, right, *e) for e in events)
    compressed(dest / "panel_events.tsv.gz", buf.getvalue())
    write_table(dest / "panel_pairs.tsv", rows)
    pairs = Counter()
    for a, b in combinations(sources, 2):
        first, second = a["ordered_id_symbol_sha256"], b["ordered_id_symbol_sha256"]
        group = "within_project" if a["project"] == b["project"] else "between_projects"
        pairs[group] += 1
        if first == second:
            pairs["same_fingerprint"] += 1
        elif results[frozenset((first, second))]["different_ids_duplicate"] > 0:
            pairs[group + "_duplicate_positive"] += 1
    for axis in axes.values():
        control, events = analyze(axis, axis)
        assert not events and not control["shared_ids_nonunique_in_either_axis"]
    def interval(field):
        values = [x[field] for x in rows]
        return [min(values), max(values)]
    summary = dict(tables=len(sources), project_labels=len({r["project"] for r in sources}),
                   fingerprints=len(axes), different_axis_pairs=len(rows),
                   different_axis_duplicate_positive=sum(x["different_ids_duplicate"] > 0 for x in rows),
                   duplicate_rebound_keys_per_different_axis_pair=interval("different_ids_duplicate"),
                   singleton_rename_ids_per_different_axis_pair=interval("singly_occurring_renamed_ids"),
                   unmatched_singleton_rename_ids_per_different_axis_pair=interval("singleton_renames_both_keys_absent_opposite"),
                   identical_axis_controls=len(axes), identical_axis_events=0,
                   sources_pair_counts=dict(sorted(pairs.items())), ordered_subset_controls=ordered_subset_control(axes))
    (dest / "panel_summary.json").write_text(json.dumps(summary, indent=2) + "\n",encoding="utf-8")
    print(json.dumps(summary, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--from-manifest", type=Path, help="Original 13-source panel manifest, beside audit_feature_panel.py")
    parser.add_argument("--output-dir", type=Path, default=HERE)
    args = parser.parse_args()
    args.output_dir.mkdir(exist_ok=True, parents=True)
    if args.from_manifest:
        extract_manifest(args.from_manifest.resolve(), args.output_dir)
    rebuild(args.output_dir if args.from_manifest else HERE, args.output_dir)


if __name__ == "__main__":
    main()
