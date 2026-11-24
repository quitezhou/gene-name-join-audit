"""Create the November follow-up ledger without changing the earlier frozen roster."""
import csv
from pathlib import Path

HERE = Path(__file__).resolve().parent
with (HERE.parent / "cohort_followup.tsv").open(encoding="utf-8-sig", newline="") as src:
    reader = csv.DictReader(src, delimiter="\t")
    fields = reader.fieldnames + ["published_abundance_after_id_alignment"]
    rows = list(reader)

assert len(rows) == len({r["pmcid"] for r in rows}) == 21
extensions = {
    "PMC11604656": "Human pseudobulk code exposed; source Visium inputs controlled, 15 join axes not traced.",
    "PMC11519383": "All 28 source axes match in order. Actual fitted gene axes and naming order unresolved.",
    "PMC10810789": "Two actual spatial H5s, archived reference counts and Figure 1 first-row rule recovered. "
                   "All 134 clinical signature scores and 132 survival tests reproduced and ID-stable. "
                   "The final cell2location signature and paired fitted abundance are unresolved."}
for row in rows:
    row["published_abundance_after_id_alignment"] = "not_adjudicated"
    if row["evidence_folder"]:
        row["evidence_folder"] = "../" + row["evidence_folder"]
    if row["pmcid"] == "PMC10017705":
        row.update(evidence_folder="marker_resolution", november_note=(
            "Published marker fractions reproduced. Partial ID repair in 1,476 source-traced cells; "
            "at least 273 detection statuses differ. Original cell2location intersection unresolved."))
    if row["pmcid"] in extensions:
        row.update(november_deep_review="1", evidence_folder="independent_cohort_extension",
                   november_note=extensions[row["pmcid"]])
    if row["pmcid"] == "PMC11116453":
        row.update(november_deep_review="1", evidence_folder="fifth_case", november_note=(
            "Reconstructed non-MT join from actual reference input axes, 81 raw reference feature tables "
            "and 40 spatial H5s; targeted reference ID independently supported by cell-count fingerprints. "
            "Observed-count Figure 3D endpoint replayed. Cell abundance was not refitted."),
            actual_c2l_identity_adjudication="different_ids_in_reconstructed_original_input")

assert sum(r["november_deep_review"] == "1" for r in rows) == 8
assert sum(r["actual_c2l_identity_adjudication"] != "indeterminate" for r in rows) == 1
with (HERE / "cohort_followup.tsv").open("w", encoding="utf-8", newline="") as dst:
    writer = csv.DictWriter(dst, fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
