from fetch_axes import *
import re
from scipy.sparse import csr_matrix
N=2300
def produce(arg):
 name,size=arg
 with Remote(endpoint('526',name),size,chunk=256*1024) as rr,h5py.File(rr,'r') as f:
  o=f['obs'];labels=o['__categories']['patient_sample'].asstr()[:];codes=o['patient_sample'][:]
  arr=np.asarray(codes).reshape(-1)
  env='Tumour' if 'Tumour' in name else 'Background'
  field='Cell types v25' if env=='Tumour' else 'Cell types v12'
  cell_label=o[field][:]
  # Choose a contiguous block with maximal single-library coverage among retained annotated cells.
  import ast
  cfg=ast.literal_eval(ast.parse((SRC/'config.py').read_text()).body[0].value)
  clnames=o['__categories'][field].asstr()[:];banned={i for i,l in enumerate(clnames) if l in cfg[env]['remove']}
  good=(cell_label>=0)&~np.isin(cell_label,list(banned))
  best=None
  for start in range(0,len(arr)-N+1,N):
   cs=arr[start:start+N]; hist=np.bincount(cs[good[start:start+N]],minlength=len(labels))
   i=max((i for i,l in enumerate(labels) if env in l),key=lambda i:hist[i])
   row=(hist[i],-start,i)
   if best is None or row>best:best=row
  matched,start,i=best;start=-start;slice_=slice(start,start+N); print('chosen',name,labels[i],matched,start,flush=True)
  ids=o['_index'].asstr()[slice_];g=f['X'];ptr=g['indptr'][start:start+N+1];lo,hi=int(ptr[0]),int(ptr[-1])
  print('extract',name,lo,hi,hi-lo,[(a,g[a].chunks) for a in g],flush=True)
  ix=g['indices'][lo:hi];v=g['data'][lo:hi];m=csr_matrix((v,ix,ptr-lo),shape=(N,g.attrs['shape'][1]))
  mask=(arr[slice_]==i)&good[slice_]; vm=m[mask]
  genes=f['var']['_index'].asstr()[:]
  cand=['ABCF2','ABCF2-1','ATXN7','ATXN7-1','DIABLO','DIABLO-1','GGT1','POLR2J3','POLR2J3-1','SCO2','RPLP0','GAPDH','RPS18','ACTB','B2M']
  with (HERE/"discordant_model_keys.tsv").open(encoding="utf8") as src_table:
   cand.extend(r["input_key"] for r in csv.DictReader(src_table,delimiter="\t"))
  cand=list(dict.fromkeys(c for c in cand if c in genes));vals=vm[:,[list(genes).index(c) for c in cand]].toarray()
  print(env,'nnz',list(zip(cand[:15],(vals[:,:15]>0).sum(0))),flush=True)
  sample=re.sub(r'Patient (\d+) (Tumour|Background) (\d+)',lambda q:'P'+q[1]+'_'+q[2][0]+q[3],labels[i])
  np.savez_compressed(SRC/(env+'_ref_window.npz'),cell=np.asarray(ids[mask],dtype=str),counts=vals,labels=np.asarray(cand),sample=sample,window_start=start,window_size=N)
  print(env,sample,len(ids[mask]),flush=True)
if __name__=='__main__':
 with ThreadPoolExecutor(2) as pool:list(pool.map(produce,[(n,v) for n,v in inventory('526').items() if n.endswith('.h5ad')]))
