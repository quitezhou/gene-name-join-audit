"""Audit order-dependent gene-name identities across a panel of 10x tables.

The input manifest names local public archives or feature files. This script
does not infer biological equivalence merely from a symbol or an Ensembl ID;
it reports the narrower software failure mode in which the same locally
unique name is rebound to another stable ID after per-dataset uniquification.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import tarfile
from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from typing import Iterable, TextIO


@dataclass(frozen=True)
class Feature:
    order: int
    stable_id: str
    symbol: str
    feature_type: str


def _feature_member(archive: tarfile.TarFile) -> tarfile.TarInfo:
    candidates = []
    for member in archive.getmembers():
        base = Path(member.name).name.lower()
        if member.isfile() and base in {
            "features.tsv",
            "features.tsv.gz",
            "genes.tsv",
            "genes.tsv.gz",
        }:
            candidates.append(member)
    if len(candidates) != 1:
        names = [member.name for member in candidates]
        raise ValueError(f"Expected one feature table in archive, found {names}")
    return candidates[0]


def open_text(path: Path) -> tuple[TextIO, str]:
    """Return a text stream and the selected member/format description."""
    if path.name.endswith((".tar.gz", ".tgz")):
        archive = tarfile.open(path, "r:gz")
        member = _feature_member(archive)
        raw = archive.extractfile(member)
        if raw is None:
            archive.close()
            raise ValueError(f"Could not extract {member.name} from {path}")
        binary: io.BufferedIOBase
        if member.name.endswith(".gz"):
            binary = gzip.GzipFile(fileobj=raw)
        else:
            binary = raw
        text = io.TextIOWrapper(binary, encoding="utf-8")
        original_close = text.close

        def close_all() -> None:
            original_close()
            archive.close()

        text.close = close_all  # type: ignore[method-assign]
        return text, member.name
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8"), "direct gzip feature table"
    return path.open("rt", encoding="utf-8"), "direct feature table"


def read_features(path: Path) -> tuple[list[Feature], str]:
    handle, selected = open_text(path)
    features: list[Feature] = []
    try:
        for line in handle:
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 2:
                continue
            feature_type = fields[2] if len(fields) >= 3 else "Gene Expression"
            if feature_type not in {"Gene Expression", ""}:
                continue
            features.append(
                Feature(
                    order=len(features),
                    stable_id=fields[0].split(".", 1)[0],
                    symbol=fields[1],
                    feature_type=feature_type or "Gene Expression",
                )
            )
    finally:
        handle.close()
    if not features:
        raise ValueError(f"No gene-expression features parsed from {path}")
    return features, selected


def by_symbol(features: Iterable[Feature]) -> dict[str, list[Feature]]:
    groups: dict[str, list[Feature]] = defaultdict(list)
    for feature in features:
        groups[feature.symbol].append(feature)
    return groups


def make_unique(features: list[Feature], join: str) -> dict[str, Feature]:
    used = {feature.symbol for feature in features}
    seen: Counter[str] = Counter()
    result: dict[str, Feature] = {}
    for feature in features:
        seen[feature.symbol] += 1
        if seen[feature.symbol] == 1:
            name = feature.symbol
        else:
            suffix = seen[feature.symbol] - 1
            name = f"{feature.symbol}{join}{suffix}"
            while name in used:
                suffix += 1
                name = f"{feature.symbol}{join}{suffix}"
            used.add(name)
        result[name] = feature
    return result


def fingerprint(features: list[Feature]) -> str:
    digest = hashlib.sha256()
    for feature in features:
        digest.update(feature.stable_id.encode())
        digest.update(b"\t")
        digest.update(feature.symbol.encode())
        digest.update(b"\n")
    return digest.hexdigest()


def summarize(features: list[Feature], selected: str) -> dict:
    groups = by_symbol(features)
    duplicates = {symbol: rows for symbol, rows in groups.items() if len(rows) > 1}
    return {
        "selected_input": selected,
        "gene_expression_features": len(features),
        "unique_stable_ids": len({feature.stable_id for feature in features}),
        "unique_symbols": len(groups),
        "duplicated_symbols": len(duplicates),
        "rows_in_duplicate_groups": sum(len(rows) for rows in duplicates.values()),
        "ordered_id_symbol_sha256": fingerprint(features),
        "duplicate_examples": [
            {
                "symbol": symbol,
                "stable_ids": [row.stable_id for row in rows],
            }
            for symbol, rows in sorted(
                duplicates.items(), key=lambda item: (-len(item[1]), item[0])
            )[:10]
        ],
    }


def compare_pair(
    first_meta: dict,
    first: list[Feature],
    second_meta: dict,
    second: list[Feature],
    join: str,
) -> dict:
    first_groups = by_symbol(first)
    second_groups = by_symbol(second)
    first_sets = {
        symbol: {row.stable_id for row in rows} for symbol, rows in first_groups.items()
    }
    second_sets = {
        symbol: {row.stable_id for row in rows}
        for symbol, rows in second_groups.items()
    }
    first_unique = make_unique(first, join)
    second_unique = make_unique(second, join)
    shared_names = first_unique.keys() & second_unique.keys()
    mismatches = []
    for name in sorted(shared_names):
        first_row = first_unique[name]
        second_row = second_unique[name]
        if first_row.stable_id == second_row.stable_id:
            continue
        symbol = first_row.symbol
        duplicate_driven = (
            len(first_groups.get(symbol, [])) > 1
            or len(second_groups.get(symbol, [])) > 1
        )
        overlap = (
            first_row.stable_id in second_sets.get(symbol, set())
            or second_row.stable_id in first_sets.get(symbol, set())
        )
        if duplicate_driven:
            mismatches.append(
                {
                    "local_unique_name": name,
                    "symbol": symbol,
                    "first_stable_id": first_row.stable_id,
                    "second_stable_id": second_row.stable_id,
                    "cross_table_id_overlap": overlap,
                }
            )
    overlap_rows = [row for row in mismatches if row["cross_table_id_overlap"]]
    return {
        "first": first_meta["label"],
        "second": second_meta["label"],
        "first_project": first_meta["project"],
        "second_project": second_meta["project"],
        "within_project": first_meta["project"] == second_meta["project"],
        "same_declared_reference": first_meta["reference"] == second_meta["reference"],
        "same_ordered_feature_fingerprint": fingerprint(first) == fingerprint(second),
        "shared_local_unique_names": len(shared_names),
        "duplicate_driven_same_name_different_id": len(mismatches),
        "duplicate_driven_with_cross_table_id_overlap": len(overlap_rows),
        "affected_symbols": len({row["symbol"] for row in mismatches}),
        "examples": overlap_rows[:8] + [
            row for row in mismatches if not row["cross_table_id_overlap"]
        ][: max(0, 8 - len(overlap_rows[:8]))],
    }


def aggregate(
    pairwise: list[dict], datasets: list[dict], summaries: dict[str, dict]
) -> dict:
    within = [row for row in pairwise if row["within_project"]]
    cross = [row for row in pairwise if not row["within_project"]]
    positives = [
        row for row in pairwise if row["duplicate_driven_same_name_different_id"] > 0
    ]
    overlap_positives = [
        row
        for row in pairwise
        if row["duplicate_driven_with_cross_table_id_overlap"] > 0
    ]
    fingerprints = {
        label: summary["ordered_id_symbol_sha256"]
        for label, summary in summaries.items()
    }
    unique_fingerprint_pairs: dict[tuple[str, str], dict] = {}
    for row in pairwise:
        first_fingerprint = fingerprints[row["first"]]
        second_fingerprint = fingerprints[row["second"]]
        if first_fingerprint == second_fingerprint:
            continue
        key = tuple(sorted((first_fingerprint, second_fingerprint)))
        unique_fingerprint_pairs.setdefault(key, row)
    provenance_flags = [
        {
            "first": row["first"],
            "second": row["second"],
            "first_reference": summaries[row["first"]]["reference"],
            "second_reference": summaries[row["second"]]["reference"],
        }
        for row in pairwise
        if row["same_ordered_feature_fingerprint"]
        and summaries[row["first"]]["reference"]
        != summaries[row["second"]]["reference"]
    ]
    return {
        "datasets": len(datasets),
        "projects": len({row["project"] for row in datasets}),
        "unique_ordered_feature_fingerprints": len(set(fingerprints.values())),
        "unique_fingerprint_pairs": len(unique_fingerprint_pairs),
        "unique_fingerprint_positive_pairs": sum(
            row["duplicate_driven_same_name_different_id"] > 0
            for row in unique_fingerprint_pairs.values()
        ),
        "pairwise_comparisons": len(pairwise),
        "within_project_pairs": len(within),
        "within_project_positive_pairs": sum(
            row["duplicate_driven_same_name_different_id"] > 0 for row in within
        ),
        "cross_project_pairs": len(cross),
        "cross_project_positive_pairs": sum(
            row["duplicate_driven_same_name_different_id"] > 0 for row in cross
        ),
        "all_positive_pairs": len(positives),
        "pairs_with_cross_table_id_overlap": len(overlap_positives),
        "maximum_duplicate_driven_mismatches_in_one_pair": max(
            (row["duplicate_driven_same_name_different_id"] for row in pairwise),
            default=0,
        ),
        "maximum_overlap_mismatches_in_one_pair": max(
            (row["duplicate_driven_with_cross_table_id_overlap"] for row in pairwise),
            default=0,
        ),
        "same_fingerprint_different_declared_reference_pairs": provenance_flags,
    }


def write_pair_csv(path: Path, pairwise: list[dict]) -> None:
    fields = [key for key in pairwise[0] if key != "examples"]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in pairwise:
            writer.writerow({key: row[key] for key in fields})


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--pair-csv", type=Path, required=True)
    args = parser.parse_args()

    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    base = args.manifest.parent
    loaded: list[tuple[dict, list[Feature]]] = []
    summaries = {}
    for metadata in manifest["datasets"]:
        path = base / metadata["local_path"]
        features, selected = read_features(path)
        loaded.append((metadata, features))
        summaries[metadata["label"]] = {
            **{key: value for key, value in metadata.items() if key != "local_path"},
            **summarize(features, selected),
        }

    pairwise = [
        compare_pair(first_meta, first, second_meta, second, "-")
        for (first_meta, first), (second_meta, second) in combinations(loaded, 2)
    ]
    pairwise.sort(
        key=lambda row: (
            -row["duplicate_driven_with_cross_table_id_overlap"],
            -row["duplicate_driven_same_name_different_id"],
            row["first"],
            row["second"],
        )
    )
    result = {
        "analysis_date": manifest["analysis_date"],
        "scope": manifest["inclusion_rule"],
        "unit_of_analysis": (
            "Feature-table pair. Pairs are not independent prevalence units; repeated "
            "samples are included to measure within-project stability."
        ),
        "dataset_summaries": summaries,
        "aggregate": aggregate(pairwise, manifest["datasets"], summaries),
        "pairwise_results": pairwise,
        "limitations": [
            "This convenience pilot is not yet a random sample of published studies.",
            "A positive pair shows a real identity hazard if tables are integrated by local names; it does not prove that the source project performed that integration.",
            "Expression and downstream consequence require the corresponding matrices and analysis code.",
        ],
    }
    args.output.write_text(json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8")
    write_pair_csv(args.pair_csv, pairwise)
    print(json.dumps(result["aggregate"], indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
