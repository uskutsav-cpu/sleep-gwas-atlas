import runpy,csv,gzip,re,collections,json,hashlib
import numpy as np,openpyxl
# Recreates exactly the fixed prior model, its training-only scaling and k-mer vocabulary.
b=runpy.run_path('work/neuronal_mpra_kmer_baseline.py')
w=openpyxl.load_workbook('work/matoba2026/table3.xlsx',read_only=True,data_only=True)
it=iter(w['a. Regulatory activity (Veh)'].iter_rows(values_only=True));next(it);loc=collections.defaultdict(set)
for r in it:
 if len(r)>9 and r[1] and r[8] and r[9]:loc[str(r[1])].add((str(r[8]),int(r[9])))
it=iter(w['b. Allelic effect (Veh)'].iter_rows(values_only=True));next(it);effects={}
for r in it:
 if len(r)>7 and r[0] is not None and isinstance(r[7],(float,int)):effects[str(r[0])]=float(r[7])
seq=collections.defaultdict(lambda:collections.defaultdict(set))
with gzip.open('work/matoba2026/barcode_statistics.gz','rt') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  m=re.search(r'_SNP_(rs\d+)_.*_(A[12])_',r['name'])
  if m and set(r['variant'])<=set('ACGT'):seq[m[1]][m[2]].add(r['variant'])
pairs=[]
for rs,p in effects.items():
 if len(loc[rs])!=1 or next(iter(loc[rs]))[0] not in ('1','2'):continue
 if len(seq[rs]['A1'])!=1 or len(seq[rs]['A2'])!=1:continue
 a=next(iter(seq[rs]['A1']));c=next(iter(seq[rs]['A2']))
 if len(a)!=150 or len(c)!=150 or sum(x!=y for x,y in zip(a,c))!=1:continue
 pairs.append((a,c,int(p<.05)))
if len(pairs)!=108 or sum(p[2] for p in pairs)!=64:raise ValueError('unexpected benchmark population')
receipt=json.load(open('outputs/stage6_matoba_reference_alignment.json'))
if receipt['status_counts']!={'aligned_reverse_complement':108}:raise ValueError('reference audit incomplete')
comp=str.maketrans('ACGT','TGCA')
def evaluate(reverse):
 X=[];y=[]
 for a,c,t in pairs:
  if reverse:a,c=a.translate(comp)[::-1],c.translate(comp)[::-1]
  X.append(np.abs(b['kcounts'](a)-b['kcounts'](c)));y.append(t)
 X=np.clip((np.stack(X)-b['mean'])/b['sd'],-20,20);y=np.array(y)
 scores=1/(1+np.exp(-np.clip(X@b['wgt']+b['intercept'],-30,30)))
 order=np.argsort(-scores,kind='stable'); labs=y[order];ss=scores[order];ends=np.r_[np.flatnonzero(np.diff(ss)!=0),len(ss)-1];tp=np.cumsum(labs)[ends];fp=(ends+1)-tp;rec=tp/sum(labs);prec=tp/(tp+fp)
 return {'auroc':float(np.trapezoid(np.r_[0,rec],np.r_[0,fp/(len(labs)-sum(labs))])),'average_precision':float(np.sum(np.diff(np.r_[0,rec])*prec))}
r={'n':108,'n_detected':64,'prevalence':64/108,'as_released_orientation':evaluate(False),'reference_plus_orientation':evaluate(True),'test':'Matoba vehicle reported allelic tests on chromosomes 1 and 2','training':'Retallick baseline training chromosomes 3-22; fixed model','warning':'Post hoc; detected/nonsignificant FDR calls, different assay and construct length, selected tests, unknown donor and sequence overlap; not score validation'}
open('outputs/stage6_matoba_external_baseline_transfer.json','w').write(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
