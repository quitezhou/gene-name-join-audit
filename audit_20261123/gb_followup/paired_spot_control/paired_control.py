"""Replay an exact-subset join on two views of the same deposited spatial spots."""
import argparse
import csv
import gzip
import hashlib
import json
from collections import Counter
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
MECHANISM = HERE.parent / "systematic_mechanism"
CASE = HERE.parent / "fifth_case"
sys.path.insert(0, str(MECHANISM))
from audit_join import assign_keys


def write_rows(path, rows):
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, rows[0], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_axis(digest):
    with gzip.open(MECHANISM / "panel_axes" / (digest + ".tsv.gz"), "rb") as f:
        text = f.read()
    assert hashlib.sha256(text).hexdigest() == digest
    return [tuple(row.split("\t")) for row in text.decode("utf-8").splitlines()]


def mapping():
    panel = json.loads((MECHANISM / "panel_summary.json").read_text())
    definition, = panel["ordered_subset_controls"]
    full = read_axis(definition["full_fingerprint"])
    subset = read_axis(definition["subset_fingerprint"])
    full_index = {record: i for i, record in enumerate(full)}
    assert len({r[0] for r in full}) == len(full) == 36601
    retained = [full_index[record] for record in subset]
    assert len(retained) == len(set(retained)) == 17943 and retained == sorted(retained)
    full_keys = list(assign_keys([r[1] for r in full], "-"))
    join = {key: i for i, key in enumerate(full_keys)}
    assert [join[full_keys[i]] for i in retained] == retained
    local = [join[key] for key in assign_keys([r[1] for r in subset], "-")]
    incorrect = [(i, j) for i, j in zip(retained, local) if i != j]
    (true_row, literal_row), = incorrect
    assert (full[true_row][0], full[literal_row][0]) == ("ENSG00000286237", "ENSG00000271147")
    return full, true_row, literal_row


def extract(cache):
    import h5py
    import numpy as np
    from scipy.sparse import csc_matrix

    full, target, literal = mapping()
    with (CASE / "input_axes.tsv").open(encoding="utf-8") as f:
        rows = [r for r in csv.DictReader(f, delimiter="\t") if r["accession_tail"] == "530"]
    output = []
    for row in rows:
        path = cache / row["file"]
        if not path.is_file():
            raise FileNotFoundError(f"{path}: cache all 40 spatial H5 files with fifth_case/replay_figure3d.py first.")
        assert hashlib.sha256(path.read_bytes()).hexdigest() == row["sha256"]
        with h5py.File(path) as f:
            m = f["matrix"]
            assert list(zip(m["features/id"].asstr()[:], m["features/name"].asstr()[:])) == full
            barcodes = m["barcodes"].asstr()[:]
            assert len(set(barcodes)) == len(barcodes) == int(row["spots"])
            matrix = csc_matrix((m["data"][:], m["indices"][:], m["indptr"][:]), shape=tuple(m["shape"][:]))
        a, b = matrix[[target, literal]].toarray()
        assert len(a) == len(b) == len(barcodes)
        # Independent source representation: sum over spots without selecting rows first.
        assert tuple(matrix.sum(axis=1).A1[[target,literal]]) == (a.sum(), b.sum())
        hist = Counter(zip(a.tolist(), b.tolist()))
        for (x, y), count in sorted(hist.items()):
            output.append(dict(section=row["file"].split("-filtered")[0], retained_id_count=x,
                               literal_full_row_count=y, spots=count))
    assert len(rows) == 40
    write_rows(HERE / "paired_histograms.tsv", output)


def summarize(rows):
    total = Counter(spots=0, neither=0, both=0, retained_only=0, literal_only=0, unequal_counts=0,
                    retained_umi=0, literal_umi=0)
    for r in rows:
        x, y, n = (int(r[k]) for k in ("retained_id_count", "literal_full_row_count", "spots"))
        assert min(x, y) >= 0 and n > 0
        category = ("both" if y else "retained_only") if x else ("literal_only" if y else "neither")
        total[category] += n
        total["spots"] += n
        total["retained_umi"] += n*x
        total["literal_umi"] += n*y
        total["unequal_counts"] += n*(x != y)
    total["either_positive"] = total["spots"] - total["neither"]
    total["discordant_detection"] = total["retained_only"] + total["literal_only"]
    assert total["either_positive"] == total["both"] + total["discordant_detection"]
    return dict(total)


def replay(path):
    import numpy as np
    _, target, literal = mapping()
    with (HERE / "paired_histograms.tsv").open(encoding="utf-8") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    sections = sorted({r["section"] for r in rows})
    assert len(sections) == 40
    by_section = [dict(section=s, **summarize([r for r in rows if r["section"] == s])) for s in sections]
    with (CASE / "input_axes.tsv").open(encoding="utf-8") as f:
        archived = {r["file"].split("-filtered")[0]:int(r["spots"]) for r in csv.DictReader(f, delimiter="\t")
                    if r["accession_tail"] == "530"}
    assert {r["section"]:r["spots"] for r in by_section} == archived
    for r in by_section:
        with np.load(CASE / "source" / (r["section"]+"-filtered_feature_bc_matrix.h5.npz")) as independent:
            assert independent["total"][[target,literal]].tolist() == [r["retained_umi"],r["literal_umi"]]
            assert independent["detected"][[target,literal]].tolist() == [
                r["both"]+r["retained_only"],r["both"]+r["literal_only"]]
    path.mkdir(parents=True, exist_ok=True)
    write_rows(path / "per_section.tsv", by_section)
    summary = dict(sections=len(sections), total=summarize(rows),
                   sections_with_detection_discordance=sum(r["discordant_detection"] > 0 for r in by_section),
                   expected_same_record_count_discrepancies=0, correctly_aligned_other_panel_rows=17942)
    (path / "summary.json").write_text(json.dumps(summary, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(summary))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--extract", type=Path, help="Read the 40 checksummed H5s and rebuild paired histograms.")
    parser.add_argument("--output-dir", type=Path, default=HERE)
    args = parser.parse_args()
    if args.extract:
        extract(args.extract)
    replay(args.output_dir)
