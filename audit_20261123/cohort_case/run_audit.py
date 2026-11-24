"""Reproduce the deposit audit. Large artifacts are cached in the OS temporary directory."""
import argparse
import os
import csv
from pathlib import Path
import subprocess
import sys
import tempfile
import zipfile

from cache_entry import cache
from derive_features import OUT
from trace_atlas import OBJECTS

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--rscript", default="Rscript", help="path to Rscript with SeuratObject/Seurat/Matrix installed")
    parser.add_argument("--refresh-inventory",action="store_true",help="re-download inventories and public feature files")
    args=parser.parse_args()
    os.environ["HWAT_CACHE"]=tempfile.gettempdir()
    def run(*parts):subprocess.run(parts,cwd=OUT,check=True)
    if args.refresh_inventory:run(sys.executable, "fetch_inputs.py")
    with (OUT/"deposit_inventory.tsv").open(encoding="utf8") as f:
        records={r["path"]:r for r in csv.DictReader(f,delimiter="\t")}
    run(sys.executable,"derive_features.py")
    run(sys.executable,"compare_panels.py")
    for name in OBJECTS:cache(records[f"Source_Data_complete/Source Data II_Seurat Files/{name}.rds"])
    run(args.rscript,"extract_objects.R",*OBJECTS)
    run(sys.executable,"trace_atlas.py")
    with zipfile.ZipFile(cache(records["se-object.visium_baseline.rds.zip"])) as archive:
        archive.extract("se-object.visium_baseline.rds",tempfile.gettempdir())
    run(args.rscript,"extract_spatial.R")
    run(sys.executable,"trace_spatial.py")
    run(sys.executable,"adjudicate_loci.py")

if __name__=="__main__":
    main()
