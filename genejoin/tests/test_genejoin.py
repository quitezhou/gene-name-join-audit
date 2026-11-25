import csv
import gzip
import subprocess
import sys
from pathlib import Path

import pytest

import genejoin
from genejoin import Axis, assign_keys, audit, recover
from genejoin.cli import main

DATA = Path(__file__).parent / "data"


def axis(rows):
    return Axis([(i, s, None) for i, s in rows])


# --- name assignment -------------------------------------------------------

def test_assign_keys_matches_anndata_rules():
    # First occurrence keeps the name; later ones skip names already present.
    assert list(assign_keys(["A", "A", "A-1", "B", "A"])) == ["A", "A-2", "A-1", "B", "A-3"]


# --- the two failure modes -------------------------------------------------

SOURCE = [("ENSG1", "X"), ("ENSG2", "X"), ("ENSG3", "Y")]


def test_filtering_before_naming_is_detected():
    a = axis(SOURCE)              # keeps both X records: X -> ENSG1, X-1 -> ENSG2
    b = axis(SOURCE[1:])          # keeps only the second: X -> ENSG2
    r = audit(a, b)
    assert not r.ok
    assert r.same_name_different_id == 1 and r.summary["different_ids_duplicate"] == 1
    assert r.same_id_different_name == 1
    ev = [e for e in r.events if e["event"] == "same_key_different_id"][0]
    assert (ev["left_id"], ev["right_id"], ev["left_key"]) == ("ENSG1", "ENSG2", "X")


def test_naming_before_filtering_is_clean():
    named = genejoin.io.name_rows(SOURCE)
    a = Axis(named)
    b = Axis(named[1:])           # inherit keys assigned on the full table
    assert audit(a, b).ok


def test_release_rename_is_detected():
    r = audit(axis([("ENSG1", "OLD"), ("ENSG9", "Z")]), axis([("ENSG1", "NEW"), ("ENSG9", "Z")]))
    assert r.same_name_different_id == 0
    assert r.same_id_different_name == 1 and r.summary["singly_occurring_renamed_ids"] == 1
    assert r.events[0]["class"] == "singleton_rename"


def test_identical_tables_are_clean():
    assert audit(axis(SOURCE), axis(SOURCE)).ok


def test_source_groups_reclassify_hidden_duplicate():
    # After filtering, B shows X only once; the full source reveals the family.
    a, b = axis([("ENSG1", "X")]), axis([("ENSG2", "X")])
    plain = audit(a, b)
    assert plain.summary["different_ids_singleton"] == 1
    b.use_source_groups(axis(SOURCE))
    a.use_source_groups(axis(SOURCE))
    assert audit(a, b).summary["different_ids_duplicate"] == 1


# --- recovery ----------------------------------------------------------------

def test_recover_uses_order():
    source = [("ENSG1", "X"), ("ENSG3", "Y"), ("ENSG2", "X"), ("ENSG4", "Z")]
    keys = ["Y", "X", "Z"]                 # bare X after Y: only ENSG2 fits in order
    assert recover(source, keys)["X"] == ("ENSG1", "ENSG2")
    assert recover(source, keys, preserve_order=True)["X"] == ("ENSG2",)


def test_recover_refuses_impossible_order():
    with pytest.raises(ValueError, match="Order check failed"):
        recover([("E1", "A"), ("E2", "B")], ["B", "A"], preserve_order=True)


def test_recover_keeps_literal_suffix_names():
    # A real symbol ending in -1 is not mistaken for a generated suffix.
    assert recover([("E1", "HLA-1"), ("E2", "HLA")], ["HLA-1"])["HLA-1"] == ("E1",)


# --- AnnData-like objects ----------------------------------------------------

class FakeAnnData:
    def __init__(self, keys, ids, symbols=None):
        pd = pytest.importorskip("pandas")
        cols = {"gene_ids": ids}
        if symbols is not None:
            cols["gene_symbols"] = symbols
        self.var = pd.DataFrame(cols, index=keys)
        self.var_names = self.var.index


def test_check_anndata():
    a = FakeAnnData(["X", "X-1", "Y"], ["ENSG1", "ENSG2", "ENSG3"], ["X", "X", "Y"])
    b = FakeAnnData(["X", "Y"], ["ENSG2", "ENSG3"], ["X", "Y"])
    r = genejoin.check_anndata(a, b)
    assert r.same_name_different_id == 1 and not r.ok


def test_check_anndata_ensembl_var_names_without_id_column():
    pd = pytest.importorskip("pandas")

    class Obj:
        var = pd.DataFrame(index=["ENSG1", "ENSG2"])
        var_names = var.index
    assert genejoin.check_anndata(Obj(), Obj()).ok


def test_check_anndata_without_ids_points_to_recover():
    a = FakeAnnData(["X"], ["ENSG1"])
    a.var = a.var.drop(columns="gene_ids")
    with pytest.raises(ValueError, match="recover"):
        genejoin.check_anndata(a, a)

@pytest.mark.parametrize("missing", [None, float("nan"), "", "  "])
def test_missing_id_is_rejected(missing):
    a = FakeAnnData(["X"], [missing])
    with pytest.raises(ValueError, match="nonempty text"):
        genejoin.check_anndata(a, a)


def test_id_fallback_validates_entire_index():
    keys = [f"ENSG{i}" for i in range(50)] + ["ENSG42-1"]
    a = FakeAnnData(keys, keys)
    a.var = a.var.drop(columns="gene_ids")
    with pytest.raises(ValueError, match="recover"):
        genejoin.check_anndata(a, a)


def test_missing_symbols_leave_mechanism_open_until_source_supplied():
    a = genejoin.axis_from_anndata(FakeAnnData(["X", "X-1"], ["E1", "E2"]))
    b = genejoin.axis_from_anndata(FakeAnnData(["X"], ["E2"]))
    result = audit(a, b)
    assert result.same_name_different_id == result.same_id_different_name == 1
    assert result.summary["different_ids_unclassified"] == 1
    assert all(e["class"] == "unknown_symbol_provenance" for e in result.events)
    for side in (a, b):
        side.use_source_groups(axis([("E1", "X"), ("E2", "X")]))
    assert audit(a, b).summary["different_ids_duplicate"] == 1


@pytest.mark.parametrize("column", ["feature_types", "feature_type"])
def test_anndata_mixed_types_keep_observed_keys_when_filtering(column):
    a = FakeAnnData(["X", "X-1"], ["Ab1", "E2"])
    a.var[column] = ["Antibody Capture", "Gene Expression"]
    with pytest.raises(ValueError, match="Non-gene"):
        genejoin.check_anndata(a, a)
    compared = genejoin.axis_from_anndata(a, mixed_modalities="genes-only")
    assert list(compared.by_key) == ["X-1"]


def test_cellxgene_gene_classes_in_backed_h5ad(tmp_path):
    ad = pytest.importorskip("anndata")
    pd = pytest.importorskip("pandas")
    ids = ["ENSG00000111640", "ENSG00000251562"]
    var = pd.DataFrame({"feature_biotype": ["gene", "gene"],
                        "feature_type": ["protein_coding", "lncRNA"],
                        "feature_name": ["GAPDH", "MALAT1"]}, index=ids)
    path = tmp_path / "cellxgene.h5ad"
    ad.AnnData(shape=(1, 2), var=var).write_h5ad(path)
    report = genejoin.check(path, path)
    assert report.ok
    assert report.summary["shared_keys"] == report.summary["shared_ids"] == 2
    assert [r.symbol for r in genejoin.load_axis(path).by_key.values()] == ["GAPDH", "MALAT1"]


def test_cellxgene_mixed_biotypes_keep_observed_keys():
    a = FakeAnnData(["X", "X-1"], ["ERCC-00002", "ENSG2"])
    a.var["feature_type"] = ["synthetic", "protein_coding"]
    a.var["feature_biotype"] = ["spike-in", "gene"]
    with pytest.raises(ValueError, match="Non-gene"):
        genejoin.check_anndata(a, a)
    compared = genejoin.axis_from_anndata(a, mixed_modalities="genes-only")
    assert [(r.id, r.key) for r in compared.by_key.values()] == [("ENSG2", "X-1")]


def test_cellxgene_and_10x_annotations_respect_both_filters():
    a = FakeAnnData(["spike", "antibody", "gene"], ["ERCC-00002", "Ab1", "ENSG2"])
    a.var["feature_biotype"] = ["spike-in", "gene", "gene"]
    a.var["feature_types"] = ["Gene Expression", "Antibody Capture", "Gene Expression"]
    assert list(genejoin.axis_from_anndata(a, mixed_modalities="genes-only").by_key) == ["gene"]


def test_partial_cellxgene_labels_report_input_error():
    a = FakeAnnData(["X", "Y"], ["E1", "E2"])
    a.var["feature_biotype"] = ["gene", None]
    with pytest.raises(ValueError, match="feature_biotype"):
        genejoin.axis_from_anndata(a, mixed_modalities="genes-only")


def test_all_spike_in_table_reports_no_genes():
    a = FakeAnnData(["spike"], ["ERCC-00002"])
    a.var["feature_biotype"] = ["spike-in"]
    with pytest.raises(ValueError, match="No gene records"):
        genejoin.axis_from_anndata(a, mixed_modalities="genes-only")


# --- command line ------------------------------------------------------------

def write_features(path, rows):
    with gzip.open(path, "wt", newline="") as f:
        csv.writer(f, delimiter="\t", lineterminator="\n").writerows([(i, s, "Gene Expression") for i, s in rows])


def test_cli_exit_codes_and_events(tmp_path, capsys):
    a, b, ev = tmp_path / "a.tsv.gz", tmp_path / "b.tsv.gz", tmp_path / "ev.tsv"
    write_features(a, SOURCE)
    write_features(b, SOURCE[1:])
    assert main(["check", str(a), str(b), "--events", str(ev)]) == 1
    assert "DISCORDANT" in capsys.readouterr().out
    rows = list(csv.DictReader(ev.open(), delimiter="\t"))
    assert {r["event"] for r in rows} == {"same_key_different_id", "same_id_different_keys"}
    assert main(["check", str(a), str(a)]) == 0
    assert main(["check", str(a), str(tmp_path / "missing.tsv")]) == 2


def test_cli_refuses_mixed_modalities(tmp_path):
    p = tmp_path / "f.tsv"
    p.write_text("ENSG1\tX\tGene Expression\nCD3\tCD3_TotalSeqB\tAntibody Capture\n")
    assert main(["check", str(p), str(p)]) == 2
    assert main(["check", str(p), str(p), "--mixed-modalities", "genes-only"]) == 0

def test_directory_features_cannot_be_replaced_by_events(tmp_path):
    source = tmp_path / "features.tsv"
    source.write_text("E1\tX\n")
    assert main(["check", str(tmp_path), str(tmp_path), "--events", str(source)]) == 2
    assert source.read_text() == "E1\tX\n"


def test_missing_optional_reader_is_an_input_error(tmp_path, monkeypatch):
    monkeypatch.setitem(sys.modules, "anndata", None)
    path = str(tmp_path / "example.h5ad")
    assert main(["check", path, path]) == 2


def test_cli_10x_folder_and_module_entry(tmp_path):
    d = tmp_path / "outs" / "filtered_feature_bc_matrix"
    d.mkdir(parents=True)
    write_features(d / "features.tsv.gz", SOURCE)
    out = subprocess.run([sys.executable, "-m", "genejoin", "check", str(tmp_path), str(tmp_path), "--json"],
                         capture_output=True, text=True)
    assert out.returncode == 0 and '"ok": true' in out.stdout


def test_cli_recover(tmp_path):
    src, keys, out = tmp_path / "s.tsv", tmp_path / "k.txt", tmp_path / "o.tsv"
    src.write_text("ENSG1\tX\nENSG3\tY\nENSG2\tX\nENSG4\tZ\n")
    keys.write_text("Y\nX\nZ\n")
    assert main(["recover", str(src), str(keys), str(out), "--order-preserved"]) == 0
    rows = list(csv.DictReader(out.open(), delimiter="\t"))
    assert [(r["final_key"], r["source_id"]) for r in rows if r["final_key"] == "X"] == [("X", "ENSG2")]


# --- the manuscript's probe-panel example (real 10x tables) -------------------

def _find_panel_tables():
    """Bundled copies in tests/data, or the manuscript's panel_axes folder via
    GENEJOIN_PANEL_AXES (file names identify the ordered source axes)."""
    ref, panel = DATA / "reference_36601.tsv.gz", DATA / "ffpe_panel_17943.tsv.gz"
    if ref.exists() and panel.exists():
        return ref, panel
    import os
    folder = os.environ.get("GENEJOIN_PANEL_AXES")
    if folder:
        root = Path(folder)
        return (root / "dd4b9d9c307903088f3f9b5231a5f72e5a79d8a02ecb3557bb6b9434597dfda2.tsv.gz",
                root / "1400d94ab85cf3feaf81e8489e11b5cde0885f1140d24ede817cf406a95e80b6.tsv.gz")
    return None, None


REF, PANEL = _find_panel_tables()
real = pytest.mark.skipif(REF is None or PANEL is None, reason="manuscript tables not available")


@real
def test_probe_panel_reassigns_armcx5_gprasp2():
    r = genejoin.check(REF, PANEL)
    assert r.same_name_different_id == 1
    ev = [e for e in r.events if e["event"] == "same_key_different_id"][0]
    assert ev["left_key"] == "ARMCX5-GPRASP2"
    assert {ev["left_id"], ev["right_id"]} == {"ENSG00000271147", "ENSG00000286237"}


@real
def test_probe_panel_recovery_with_order():
    source = genejoin.read_feature_table(REF)
    panel = genejoin.read_feature_table(PANEL)
    keys = [k for _, _, k in genejoin.io.name_rows(panel)]
    plain = recover(source, keys)
    ordered = recover(source, keys, preserve_order=True)
    assert sum(len(v) > 1 for v in plain.values()) == 9
    ambiguous = [k for k, v in ordered.items() if len(v) > 1]
    assert ambiguous == ["MATR3"]
    assert ordered["ARMCX5-GPRASP2"] == ("ENSG00000286237",)
