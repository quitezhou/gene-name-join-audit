"""Extract the prespecified whole published Tables S10 and S11 from the publisher PDF."""
import argparse
import csv
import hashlib
from pathlib import Path
import re

from pypdf import PdfReader

HERE = Path(__file__).resolve().parent


def main(pdf):
    assert hashlib.sha256(pdf.read_bytes()).hexdigest() == "d8a32487c1b1608a5d6dd0531030fdd91bcd18c2677db0317a1e04bf4dd7726b"
    reader = PdfReader(pdf)
    page = reader.pages[26].extract_text()
    genes = page[page.index("TMEM176A"):page.index("CXCL-\nifibro")].split()
    assert len(set(genes)) == len(genes) == 65 and genes[-1] == "SELENOP"
    (HERE / "kidney_S10_symbols.txt").write_text("\n".join(genes)+"\n", encoding="utf8")
    page = reader.pages[27].extract_text(extraction_mode="layout")
    records = []
    for line in page.splitlines():
        if not re.match(r"^\s*\d+\.\d+", line):
            continue
        m = re.fullmatch(r"\s*(\d+\.\d+)\s+(\d+)\s+(\d+)\s+(0\.\d+)\s+(\d+\.\d+)\s+(\d+)\s+(\d+)\s*(0\.\d+)\s*", line)
        if m is None:
            raise ValueError(line)
        for start in [1,5]:
            rec = [m.group(i) for i in range(start,start+4)]
            assert int(rec[1])+int(rec[2]) == 134
            records.append(rec)
    assert sorted(int(r[1]) for r in records) == list(range(1,133))
    with (HERE / "kidney_S11_published.tsv").open("w", encoding="utf8",newline="") as f:
        w = csv.writer(f,delimiter="\t",lineterminator="\n")
        w.writerow(["threshold","n_low","n_high","p"])
        w.writerows(sorted(records,key=lambda r:int(r[1])))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("pdf",type=Path)
    main(ap.parse_args().pdf)
