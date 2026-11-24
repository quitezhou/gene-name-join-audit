"""Apply the standalone keyed-input checker to the archived lung modeling axes."""
import argparse
import csv
import gzip
import hashlib
import json
import subprocess
import sys
import tempfile
from pathlib import Path

from audit_join import Axis, parse

HERE = Path(__file__).resolve().parent
CASE = HERE.parent / "fifth_case"


def records(path):
    with path.open(encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def source_axis(manifest, accession_tail, libraries):
    rows = [r for r in manifest if r["accession_tail"] == accession_tail]
    axes = {r["axis_id"] for r in rows}
    assert len(rows) == libraries and len(axes) == 1
    digest = axes.pop()
    with gzip.open(CASE / "source" / f"{digest}.tsv.gz", "rt", encoding="utf-8", newline="") as f:
        entries = parse(f)
    assert all(len(entries) == int(r["features"]) for r in rows)
    return entries, CASE / "source" / f"{digest}.tsv.gz"


def keyed_export(path, rows):
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(["stable_id", "symbol", "final_key"])
        writer.writerows((r.id, r.symbol, r.key) for r in rows)


def main(output):
    manifest = records(CASE / "input_axes.tsv")
    sc_entries, sc_src = source_axis(manifest, "526", 81)
    sc = Axis(sc_entries)
    spatial, st_src = source_axis(manifest, "530", 40)
    non_mt = Axis([r for r in spatial if not r[1].startswith("MT-")])
    # Reference keys were assigned on the full raw table before subsetting.
    keys = [[r["_index"] for r in records(p)] for p in sorted((CASE / "source").glob("*.h5ad.var.tsv"))]
    assert len(keys) == 2 and keys[0] == keys[1] and len(set(keys[0])) == len(keys[0])
    ref = [sc.by_key[k] for k in keys[0]]
    with tempfile.TemporaryDirectory() as temp:
        folder = Path(temp)
        sc_tsv, st_tsv, events = (folder / name for name in ("sc.tsv", "spatial.tsv", "events.tsv"))
        keyed_export(sc_tsv, ref)
        keyed_export(st_tsv, non_mt.by_key.values())
        cmd = [sys.executable, str(HERE / "audit_join.py"), str(sc_tsv), str(st_tsv),
               "--left-keyed", "--right-keyed", "--events", str(events)]
        def replay(extra):
            run = subprocess.run(cmd + extra, check=True, capture_output=True, text=True)
            return json.loads(run.stdout), records(events)
        final_only, final_events = replay([])
        summary, rows = replay(["--left-source-features", str(sc_src), "--right-source-features", str(st_src)])
    kind = {"duplicate_associated": "duplicate_associated", "singleton_on_both": "singleton_name"}
    found = {r["left_key"]: (r["left_symbol"], r["right_symbol"], r["left_id"], r["right_id"], kind[r["class"]])
             for r in rows if r["event"] == "same_key_different_id"}
    original = {r["input_key"]: (r["sc_base"], r["vis_base"], r["sc_source_id"], r["vis_id"], r["cause"])
                for r in records(CASE / "discordant_model_keys.tsv")}
    assert len(found) == len(original) == 221 and found == original
    final_mismatch = {r["left_key"]: r for r in final_events if r["event"] == "same_key_different_id"}
    assert set(final_mismatch) == set(found)
    for k, (sc_base, vis_base, sc_id, vis_id, _) in original.items():
        r = final_mismatch[k]
        assert (r["left_symbol"], r["right_symbol"], r["left_id"], r["right_id"]) == (sc_base, vis_base, sc_id, vis_id)
    reclassified = [k for k in sorted(found) if kind[final_mismatch[k]["class"]] != found[k][-1]]
    assert reclassified == ["SCO2"] and (final_only["different_ids_duplicate"],final_only["different_ids_singleton"]) == (5,216)
    assert (summary["shared_keys"], summary["different_ids_duplicate"], summary["different_ids_singleton"],
            summary["shared_ids"], summary["shared_id_different_keys"],
            summary["shared_ids_nonunique_in_either_axis"]) == (24769, 6, 215, 25302, 754, 0)
    assert (len(ref), len(non_mt.by_key)) == (25718, 36588)
    canonical = "".join("\t".join((k, *v))+"\n" for k, v in sorted(found.items()))
    report = dict(study="De Zuani et al., doi:10.1038/s41467-024-48700-8",
                  comparison="Independent implementation on the same reconstructed original input axes",
                  shared_reference_input_keys=len(ref), non_mitochondrial_spatial_keys=len(non_mt.by_key),
                  shared_keys=summary["shared_keys"], same_key_different_ids=len(found),
                  duplicate_associated=summary["different_ids_duplicate"], singleton_name=summary["different_ids_singleton"],
                  shared_ids=summary["shared_ids"], shared_ids_missed=summary["shared_id_different_keys"],
                  archived_mismatch_rows_exactly_matched=len(found),
                  final_axes_without_source_groups=dict(duplicate_associated=final_only["different_ids_duplicate"],
                      singleton_on_final_axes=final_only["different_ids_singleton"], keys_identified_as_source_duplicates_only_with_provenance=reclassified),
                  comparison_fields=["input_key", "sc_base", "vis_base", "sc_source_id", "vis_id", "cause"],
                  sorted_row_digest_sha256=hashlib.sha256(canonical.encode()).hexdigest())
    output.write_text(json.dumps(report, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    arg = argparse.ArgumentParser(description=__doc__)
    arg.add_argument("--output", type=Path, default=HERE / "lung_join_check.json")
    main(arg.parse_args().output)
