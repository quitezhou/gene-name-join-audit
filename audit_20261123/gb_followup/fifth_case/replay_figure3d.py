"""Replay the published colocalisation endpoint on source H5 counts, before/after re-keying."""
import csv
import hashlib
import json
import tempfile
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csc_matrix
from compare_axes import read_axis, key_map
from fetch_axes import SRC, HERE, endpoint, get_all

CACHE = Path(tempfile.gettempdir()) / 'nsclc_ref_check'
CACHE.mkdir(exist_ok=True)

def cached_input(row):
    path = CACHE / row.file
    if not path.exists() or hashlib.sha256(path.read_bytes()).hexdigest() != row.sha256:
        content = get_all(endpoint('530', row.file), row.bytes)
        assert hashlib.sha256(content).hexdigest() == row.sha256
        path.write_bytes(content)
    return path

def data_from_slide(job):
    row, masks, gene_rows, gene_order, pairs = job
    path = cached_input(row)
    with h5py.File(path) as h:
        g = h['matrix']
        m = csc_matrix((g['data'][:],g['indices'][:],g['indptr'][:]),shape=tuple(g['shape'][:])).tocsr()
        barcodes=g['barcodes'].asstr()[:]; ids=g['features']['id'].asstr()[:]
    gene_counts = m[gene_rows,:].toarray()
    sample = row.file.split('-filtered')[0]
    sample_rows = []
    spot_qc = {}; qc_totals = {}
    for method, flag in masks.items():
        total = np.asarray(m[flag,:].sum(axis=0)).ravel()
        keep = total >= 800
        spot_qc[method] = keep; qc_totals[method] = total
        x = gene_counts[:,keep]
        active = x > np.median(x,axis=1,keepdims=True)
        on = {gene: active[i,:] for i,gene in enumerate(gene_order)}
        for pair in pairs:
            l,r = pair.split('_')
            sample_rows.append([sample, method, pair, int((on[l]&on[r]).sum()),int(keep.sum())])
    qcs = [[sample,method,int(keep.sum()),int((keep & ~spot_qc['source']).sum()),int((~keep & spot_qc['source']).sum())] for method,keep in spot_qc.items()]
    strict='six_duplicate_groups_by_id'; changes=[]
    altered=np.flatnonzero(masks['source']!=masks[strict])
    for spot in np.flatnonzero(spot_qc['source']!=spot_qc[strict]):
        gene_counts=m[altered,spot].toarray().ravel()
        for j,gene_row in enumerate(altered):
            changes.append([sample,barcodes[spot],int(qc_totals['source'][spot]),int(qc_totals[strict][spot]),ids[gene_row],bool(masks[strict][gene_row]),int(gene_counts[j])])
    return sample_rows,qcs,changes

def main():
    manifest = pd.read_csv(HERE/'input_axes.tsv',sep='\t',dtype={'accession_tail':str})
    sc = read_axis(manifest.loc[manifest.accession_tail=='526','axis_id'].iloc[0]); vis = read_axis(manifest.loc[manifest.accession_tail=='530','axis_id'].iloc[0]); visit = manifest[manifest.accession_tail=='530']
    sc = key_map(sc)
    allowed = pd.read_csv(next(SRC.glob('*Tumour*.h5ad.var.tsv')),sep='\t',keep_default_na=False).iloc[:,0]
    sc = sc.loc[allowed]
    vis['row'] = np.arange(len(vis))
    spatial = key_map(vis[~vis.symbol.str.startswith('MT-')]); ref_ids = set(sc.id)
    original = spatial.index.isin(sc.index); exact = spatial.id.isin(ref_ids).to_numpy()
    bad = pd.read_csv(HERE/'discordant_model_keys.tsv',sep='\t'); six=set(bad.loc[bad.cause=='duplicate_associated','sc_base']); assert len(six)==6
    local = np.where(spatial.symbol.isin(six),exact,original)
    def all_rows(keep):
        mask=np.zeros(len(vis),dtype=bool); mask[spatial.row.to_numpy()[keep]]=True;return mask
    masks={name:all_rows(keep) for name,keep in [('source',original),('six_duplicate_groups_by_id',local),('all_shared_ids',exact)]}
    pub = pd.read_csv(SRC/'LRcolocalisation_SpotPercentage_TvsB.csv',index_col=0,header=[0,1])
    pairs = list(pub.index)
    go = sorted({g for p in pairs for g in p.split('_')}); gene_rows = spatial.loc[go,'row'].to_list()
    assert all(spatial.loc[g,'id']==sc.loc[g,'id'] for g in go)
    result=[]; quality=[]; traces=[]
    task=[(r,masks,gene_rows,go,pairs) for r in visit.itertuples(index=False)]
    with ThreadPoolExecutor(4) as ex:
        for sample_rows,qcs,changes in ex.map(data_from_slide,task):result.extend(sample_rows); quality.extend(qcs); traces.extend(changes)
    pd.DataFrame(traces,columns=['sample','barcode','original_UMI','corrected_UMI','changed_membership_ID','present_after_fix','this_gene_UMI']).to_csv(HERE/'six_group_changed_spots.tsv',sep='\t',index=False)
    qc=pd.DataFrame(quality,columns=['sample','method','passing_spots','gained_vs_original','lost_vs_original']);qc.to_csv(HERE/'figure3d_qc.tsv',index=False,sep='\t')
    d=pd.DataFrame(result,columns=['sample','method','pair','both','passing'])
    sdrf=pd.read_csv(SRC/'E-MTAB-13530.sdrf.txt',sep='\t'); sdrf['id']=sdrf['Scan Name'].str.extract(r'LNGsp(\d+)',expand=False)
    md=pd.read_csv(SRC/'Spatial_metadata.csv');md['id']=md.CellRanger.str.extract(r'LNGsp(\d+)',expand=False)
    key=md.merge(sdrf[['id','Source Name']].drop_duplicates(),on='id',validate='one_to_one');assert len(key)==40
    d=d.merge(key,left_on='sample',right_on='Source Name',validate='many_to_one');d=d[d.Tissue!='Healthy']
    gp=d.groupby(['Patient_number','Tissue','pair','method'])[['both','passing']].sum()
    gp.to_csv(HERE/'figure3d_spot_counts.tsv',sep='\t')
    obs=(100 * gp.both / gp.passing).unstack('method'); obs['author_pct']=[pub.loc[p,(pt,tis)] for pt,tis,p in obs.index];obs['source_error_pp']=obs.source-obs.author_pct
    for name in ['six_duplicate_groups_by_id','all_shared_ids']:obs[name+'_delta_pp']=obs[name]-obs.source
    obs.to_csv(HERE/'figure3d_endpoint.tsv',sep='\t');print('Counts, original residual, six/full largest pp delta:',len(obs),[float(obs[c].abs().max()) for c in ['source_error_pp','six_duplicate_groups_by_id_delta_pp','all_shared_ids_delta_pp']]);print(qc.groupby('method').sum(numeric_only=True));print('variable-axis sizes', {k:int(v.sum()) for k,v in masks.items()})
if __name__=='__main__': main()
