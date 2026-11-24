"""Install the archived wheel in a temporary directory and audit the real probe panel."""

import csv
import gzip
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import runpy


def check(root):
    study = root / "audit_20261123" / "gb_followup"
    evidence = study / "systematic_mechanism"
    control, = json.loads((evidence / "panel_summary.json").read_text())["ordered_subset_controls"]
    axes = evidence / "panel_axes"
    full = axes / (control["full_fingerprint"] + ".tsv.gz")
    subset = axes / (control["subset_fingerprint"] + ".tsv.gz")
    wheel, = (root / "genejoin" / "dist").glob("*.whl")
    with tempfile.TemporaryDirectory(prefix="genejoin_install_") as tmp:
        work = Path(tmp)
        target = work / "site"
        subprocess.run(
            [sys.executable, "-m", "pip", "install", "--no-index", "--no-deps", "--no-compile",
             "--disable-pip-version-check", "--target", str(target), str(wheel)],
            cwd=work, check=True, stdout=subprocess.DEVNULL,
        )
        env = dict(os.environ, PYTHONPATH=str(target) + os.pathsep + os.environ.get("PYTHONPATH", ""))
        subprocess.run([sys.executable, "-s", "-c",
                        "from pathlib import Path; import genejoin, sys; "
                        "assert Path(genejoin.__file__).resolve().is_relative_to(Path(sys.argv[1]).resolve())",
                        str(target)], cwd=work, env=env, check=True)

        def run(*args):
            return subprocess.run([sys.executable, "-s", "-m", "genejoin", *map(str, args)],
                                  cwd=work, env=env, capture_output=True, text=True)

        version = run("--version")
        assert version.returncode == 0 and version.stdout.strip() == "genejoin 0.1.0", version.stderr
        clean = run("check", subset, subset)
        assert clean.returncode == 0, clean.stderr
        events = work / "events.tsv"
        discordant = run("check", full, subset, "--events", events)
        assert discordant.returncode == 1, discordant.stderr
        with events.open(encoding="utf-8", newline="") as fh:
            observed = list(csv.DictReader(fh, delimiter="\t"))
        fields = ("event", "left_key", "right_key", "left_id", "right_id")
        def signature(records):
            return sorted(tuple(r[k] for k in fields) for r in records)
        assert signature(observed) == signature(control["independent_discrepancies"])
        assert len(observed) == 2 and {r["left_symbol"] for r in observed} == {"ARMCX5-GPRASP2"}

        with gzip.open(subset, "rt", encoding="utf-8", newline="") as fh:
            retained = [row[:2] for row in csv.reader(fh, delimiter="\t")]
        # Generate the expected labels with the original independent panel checker.
        assign_keys = runpy.run_path(str(evidence / "audit_join.py"))["assign_keys"]
        keys = list(assign_keys([r[1] for r in retained], "-"))
        assert len(set(keys)) == len(keys) == 17943
        key_file = work / "keys.txt"
        key_file.write_text("\n".join(keys) + "\n", encoding="utf-8")
        output = work / "candidates.tsv"
        recovered = run("recover", full, key_file, output, "--order-preserved")
        assert recovered.returncode == 0, recovered.stderr
        expected_ids = dict(zip(keys, [r[0] for r in retained]))
        with output.open(encoding="utf-8", newline="") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        matched = {r["final_key"] for r in rows if r["source_id"] == expected_ids[r["final_key"]]}
        ambiguous = {r["final_key"] for r in rows if int(r["candidate_count"]) > 1}
        assert matched == set(keys) and ambiguous == {"MATR3"}
        print("Installed genejoin 0.1.0 wheel reproduced both join discrepancies and the MATR3 ambiguity.")


if __name__ == "__main__":
    check(Path(__file__).resolve().parent)
