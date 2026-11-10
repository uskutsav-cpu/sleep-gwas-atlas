import io,tarfile,re,hashlib,csv,json
from collections import defaultdict
import h5py,numpy as np
from scipy.sparse import csc_matrix
P='work/sleep_sc/GSE275071_RAW.tar'; genes=['FOXP2','SMAD5','RBM6','MEIS1']; markers=['CD14','LST1','CD3D','CD3E']
acc=defaultdict(lambda: {'cells':0,'detect':dict.fromkeys(genes,0),'lanes':0})
with tarfile.open(P) as tar:
 members=tar.getmembers(); assert len(members)==32
 for member in members:
  m=re.search(r'Lane([1-4])_(HS|SR)-([1-4])_',member.name)
  if not m: raise ValueError(member.name)
  lane,group,donor=m.groups()
  with h5py.File(io.BytesIO(tar.extractfile(member).read()),'r') as h:
   g=h['matrix']; names=np.array([x.decode() for x in g['features']['name'][:]])
   ix={}
   for gene in genes+markers:
    indices=np.flatnonzero(names==gene)
    if len(indices)!=1: raise ValueError(f'{gene}: {len(indices)}')
    ix[gene]=indices[0]
   a=csc_matrix((g['data'][:],g['indices'][:],g['indptr'][:]),shape=g['shape'][:])
   ncounts=np.asarray(a.sum(axis=0)).ravel(); detected=np.diff(a.indptr)
   mito=np.asarray(a[np.char.startswith(names.astype(str),'MT-'),:].sum(axis=0)).ravel()
   qc=(ncounts>=500)&(detected>=200)&(mito/np.maximum(ncounts,1)<.2)
   marker={gene:np.asarray(a[ix[gene],:].todense()).ravel() for gene in markers}
   types={'monocyte_marker_positive':qc&(marker['CD14']>0)&(marker['LST1']>0)&(marker['CD3D']==0)&(marker['CD3E']==0),
          't_marker_positive':qc&(marker['CD3D']>0)&(marker['CD3E']>0)&(marker['CD14']==0)&(marker['LST1']==0),
          'all_qc_cells':qc,
          'cd14_marker_positive_relaxed':qc&(marker['CD14']>0),
          'cd3d_marker_positive_relaxed':qc&(marker['CD3D']>0)}
   for typ,mask in types.items():
    d=acc[(group,donor,typ)];d['lanes']+=1;d['cells']+=int(mask.sum())
    for gene in genes:d['detect'][gene]+=int(np.count_nonzero(np.asarray(a[ix[gene],:].todense()).ravel()[mask]))
assert len(acc)==40 and all(a['lanes']==4 for a in acc.values())
rows=[]
for group in ['HS','SR']:
 for typ in ['all_qc_cells','cd14_marker_positive_relaxed','cd3d_marker_positive_relaxed']:
  samples=[v for (g,d,t),v in acc.items() if g==group and t==typ]
  for gene in genes:
   fractions=[v['detect'][gene]/v['cells'] for v in samples if v['cells']]
   rows.append({'group':group,'marker_defined_subset':typ,'gene':gene,'donors':len(fractions),'total_cells':sum(v['cells'] for v in samples),'mean_donor_detection_fraction':float(np.mean(fractions)),'min_donor_fraction':min(fractions),'max_donor_fraction':max(fractions)})
with open('outputs/stage6_gse275071_gate_sensitivity.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
receipt={'source':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE275nnn/GSE275071/suppl/GSE275071_RAW.tar','sha256':hashlib.sha256(open(P,'rb').read()).hexdigest(),'h5_files':32,'donors_per_group':4,'wells_per_donor':4,'qc':'UMI >= 500, detected genes >= 200, mitochondrial fraction < 0.2','summary':'unweighted mean, minimum and maximum of donor-level detected-cell fractions','warning':'marker-gated descriptive detection only; no inferential testing, cell type validation, insomnia/frailty or genotype'}
open('outputs/stage6_gse275071_gate_sensitivity.json','w').write(json.dumps(receipt,indent=2)+'\n')
for row in rows:print(row)
