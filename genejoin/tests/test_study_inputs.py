"""Regressions on frozen public feature metadata and two independently released H5ADs.

GENEJOIN_STUDY_DIR points to the extracted audit_20261123/gb_followup directory.
The two optional H5AD variables below point to the specified downloads. The tests
do not make network requests or need expression matrices for the study replays.

GENEJOIN_PBMC3K_H5AD: https://falexwolf.de/data/pbmc3k_raw.h5ad
GENEJOIN_SYMBOL_ONLY_H5AD: https://raw.githubusercontent.com/chanzuckerberg/cellxgene/main/example-dataset/pbmc3k.h5ad
"""

import csv
import gzip
import hashlib
import json
import os
from itertools import combinations
from pathlib import Path

import pytest

import genejoin
from genejoin.cli import main
from genejoin.io import name_rows


def provided(name):
    value = os.environ.get(name)
    if not value:
        pytest.skip(f"Set {name} to run this public-data regression")
    path = Path(value)
    assert path.exists(), f"{name} points to a missing input: {path}"
    return path


def records(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_keyed(path, rows):
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(("stable_id", "symbol", "final_key"))
        writer.writerows(rows)


def test_released_lung_model_axes(tmp_path, capsys):
    case = provided("GENEJOIN_STUDY_DIR") / "fifth_case"
    manifest = records(case / "input_axes.tsv")

    def source(tail, libraries):
        matches = [row["axis_id"] for row in manifest if row["accession_tail"] == tail]
        assert len(matches) == libraries and len(set(matches)) == 1
        return case / "source" / (matches[0] + ".tsv.gz")

    sc, spatial = source("526", 81), source("530", 40)
    axes = [[row["_index"] for row in records(path)]
            for path in sorted((case / "source").glob("*.h5ad.var.tsv"))]
    assert len(axes) == 2 and axes[0] == axes[1] and len(axes[0]) == 25718

    # The reference was named on the full source before final subsetting.
    sc_rows = name_rows(genejoin.read_feature_table(sc))
    sc_by_key = {row[2]: row for row in sc_rows}
    left, right, events = (tmp_path / name for name in ("sc.tsv", "spatial.tsv", "events.tsv"))
    write_keyed(left, (sc_by_key[key] for key in axes[0]))
    spatial_rows = [r for r in genejoin.read_feature_table(spatial) if not r[1].startswith("MT-")]
    write_keyed(right, name_rows(spatial_rows))

    assert main(["check", str(left), str(right), "--left-source", str(sc), "--right-source",
                 str(spatial), "--json", "--events", str(events)]) == 1
    report = json.loads(capsys.readouterr().out)
    assert (report["same_name_different_id"], report["same_id_different_name"]) == (221, 754)
    summary = report["summary"]
    assert (summary["shared_keys"], summary["shared_ids"]) == (24769, 25302)
    assert (summary["different_ids_duplicate"], summary["different_ids_singleton"]) == (6, 215)

    classes = {"duplicate_associated": "duplicate_associated", "singleton_on_both": "singleton_name"}
    actual = {(e["left_key"], e["left_symbol"], e["right_symbol"], e["left_id"], e["right_id"],
               classes[e["class"]]) for e in records(events) if e["event"] == "same_key_different_id"}
    fields = ("input_key", "sc_base", "vis_base", "sc_source_id", "vis_id", "cause")
    expected = {tuple(row[f] for f in fields) for row in records(case / "discordant_model_keys.tsv")}
    assert len(actual) == 221 and actual == expected


def test_withheld_blood_ids_on_all_axes_and_primary_pairs():
    study = provided("GENEJOIN_STUDY_DIR")
    counts = records(study / "prospective_holdout" / "axis_counts.tsv.gz")
    samples = records(study / "prospective_holdout" / "samples.tsv")
    source = [(r["id"], r["symbol"]) for r in counts]
    reference = genejoin.Reference(source)
    bounds = study / "prejoin_bounds"
    expected_axes = {(r["sample"], int(r["threshold"])): r for r in records(bounds / "axes.tsv")
                     if r["dataset"] == "PBMC"}
    assert len(samples) == 16 and len(expected_axes) == 48
    primary = {}
    for threshold in (1, 3, 10):
        for sample in samples:
            label = sample["sample"]
            survivors = [item for item, r in zip(source, counts) if int(r[label + "_cells"]) >= threshold]
            gold = {key: fid for fid, _, key in name_rows(survivors)}
            keys = list(gold)
            unordered = reference.candidates(keys)
            ordered = reference.candidates(keys, preserve_order=True)
            assert all(gold[key] in ordered[key] and set(ordered[key]) <= set(unordered[key]) for key in keys)
            expected = expected_axes[label, threshold]
            assert len(keys) == int(expected["rows"])
            assert sum(len(ids) > 1 for ids in unordered.values()) == int(expected["ambiguous_without_order"])
            assert sum(len(ids) > 1 for ids in ordered.values()) == int(expected["ambiguous_with_order"])
            if threshold == 1:
                primary[label] = (gold, ordered)

    expected_pairs = {(r["left"], r["right"]): r for r in records(bounds / "pairs.tsv")}
    assert len(expected_pairs) == 56
    totals = []
    for condition in ("ctrl", "stim"):
        labels = [r["sample"] for r in samples if r["condition"] == condition]
        assert len(labels) == 8
        differences = certifications = 0
        for left, right in combinations(labels, 2):
            gold_a, bounds_a = primary[left]
            gold_b, bounds_b = primary[right]
            shared = gold_a.keys() & gold_b.keys()
            true = {k for k in shared if gold_a[k] != gold_b[k]}
            different = {k for k in shared if set(bounds_a[k]).isdisjoint(bounds_b[k])}
            same = {k for k in shared if len(bounds_a[k]) == 1 and bounds_a[k] == bounds_b[k]}
            expected = expected_pairs[left, right]
            assert different <= true and same.isdisjoint(true)
            assert (len(shared), len(same), len(different), len(true)) == tuple(
                int(expected[f]) for f in ("shared", "certified_same", "certified_different", "actual_differences"))
            assert "|".join(sorted(true - different)) == expected["unresolved_actual_keys"]
            differences += len(true)
            certifications += len(different)
        totals.append((differences, certifications))
    assert totals == [(51, 14), (55, 13)]


def test_public_pbmc3k_h5ad_against_independent_10x_axis(capsys):
    study = provided("GENEJOIN_STUDY_DIR")
    h5ad = provided("GENEJOIN_PBMC3K_H5AD")
    anndata = pytest.importorskip("anndata")
    assert hashlib.sha256(h5ad.read_bytes()).hexdigest() == "89a96f1beaa2dd83a687666d3f19a4513ac27a2a2d12581fcd77afed7ea653a1"
    source = study / "systematic_mechanism" / "panel_axes" / "8778dd78085086be6dcafa8334f95f7cf76b25526e8ef2b63186db5127a1492c.tsv.gz"
    archived = genejoin.load_axis(source)
    adata = anndata.read_h5ad(h5ad, backed="r")
    try:
        assert adata.shape == (2700, 32738)
        assert list(adata.var_names) == list(archived.by_key)
        assert list(adata.var["gene_ids"]) == [row.id for row in archived.by_key.values()]
    finally:
        adata.file.close()
    assert main(["check", str(h5ad), str(source), "--json"]) == 0
    report = json.loads(capsys.readouterr().out)
    assert report["summary"]["shared_ids"] == report["summary"]["shared_keys"] == 32738
    assert report["same_name_different_id"] == report["same_id_different_name"] == 0


def test_public_processed_h5ad_without_ids_fails_closed(capsys):
    h5ad = provided("GENEJOIN_SYMBOL_ONLY_H5AD")
    pytest.importorskip("anndata")
    assert hashlib.sha256(h5ad.read_bytes()).hexdigest() == "0db367b991dd95809732b218539ede489bea99113807f62ebd7ccc970025fe38"
    assert main(["check", str(h5ad), str(h5ad), "--json"]) == 2
    captured = capsys.readouterr()
    assert not captured.out and "stable ID" in captured.err
