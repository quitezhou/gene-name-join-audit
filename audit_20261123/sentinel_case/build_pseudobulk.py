"""Build pseudobulk counts from public GEO MEX, applying two numerical QC steps.

Usage: python build_pseudobulk.py --cache PATH [--workers 4]
Run in an individual writable cache; retained sample-level counts can be audited
from cached GEO files. The deposited matrices have already been filtered by 10x.
"""
import argparse
import csv
import gzip
import hashlib
import json
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from scipy.io import mmread


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SOURCES = ROOT / 'panel_sources'

def process(task):
    spec, cache = task
    filename, expected = spec['Name'], int(spec['Size'])
    acc, group, number, _ = filename.split('_')
    source = SOURCES / filename
    url = f'https://ftp.ncbi.nlm.nih.gov/geo/samples/{acc[:-3]}nnn/{acc}/suppl/{filename}'
    if not source.is_file():
        source = cache / filename
    if not source.exists() or source.stat().st_size != expected:
        partial = source.with_suffix('.part')
        subprocess.run(['curl.exe' if os.name=='nt' else 'curl', '-fLsS', '--connect-timeout', '15',
                        '--max-time', '240', '--retry', '3', url, '-o', str(partial)], check=True)
        if partial.stat().st_size != expected:
            raise ValueError(f'GEO size mismatch for {filename}')
        partial.replace(source)
    with gzip.open(source, 'rb') as f:
        mat = mmread(f).tocsr()
    if np.any(mat.data <= 0):
        raise ValueError('Unexpected explicit zero or negative count')
    genes = np.diff(mat.indptr) >= 5
    qc = mat[genes]
    spots = np.asarray(qc.sum(axis=0)).ravel() >= 500
    counts = np.zeros(mat.shape[0],dtype='int64')
    counts[genes] = np.asarray(qc[:,spots].sum(axis=1)).ravel()
    raw = np.asarray(mat.sum(axis=1)).ravel()
    metrics = dict(sample=f'{group}_{number}',accession=acc,url=url,bytes=expected,
                   sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                   raw_spots=mat.shape[1],qc_spots=int(spots.sum()),raw_UMI=int(raw.sum()),qc_UMI=int(counts.sum()))
    print(f"{acc} {group}_{number} {metrics['qc_UMI']}",flush=True)
    return (acc,counts,raw,metrics)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache',type=Path,required=True)
    parser.add_argument('--workers',type=int,default=4)
    args=parser.parse_args()
    args.cache.mkdir(parents=True,exist_ok=True)
    with (SOURCES/'GSE278936_filelist.txt').open() as f:
        specs=[r for r in csv.DictReader(f,delimiter='\t') if r['Name'].endswith('_matrix.mtx.gz')
               and r['Name'].split('_')[1] in ('TRNA','NEADT')]
    assert len(specs)==39
    rows=[]
    if (HERE/'matrix_manifest.json').is_file() and (HERE/'pseudobulk_counts.npz').is_file():
        previous=json.loads((HERE/'matrix_manifest.json').read_text())
        with np.load(HERE/'pseudobulk_counts.npz') as saved:
            for m in previous:
                acc=m['accession']
                rows.append((acc,saved[f'{acc}_qc'], saved[f'{acc}_raw'],m))
        known={r[0] for r in rows}
        specs=[s for s in specs if s['Name'].split('_')[0] not in known]
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = [(s['Name'], pool.submit(process,(s,args.cache))) for s in specs]
        unavailable = []
        for name, future in futures:
            try:
                rows.append(future.result())
            except (subprocess.CalledProcessError, OSError) as exc:
                unavailable.append(dict(filename=name, reason=str(exc)))
    (HERE/'unavailable.json').write_text(json.dumps(unavailable,indent=2)+'\n')
    np.savez_compressed(HERE/'pseudobulk_counts.npz',**{f'{acc}_{k}':v for acc,q,r,_ in rows for k,v in [('qc',q),('raw',r)]})
    (HERE/'matrix_manifest.json').write_text(json.dumps([m for *_,m in rows],indent=2)+'\n')

if __name__=='__main__':
    main()



