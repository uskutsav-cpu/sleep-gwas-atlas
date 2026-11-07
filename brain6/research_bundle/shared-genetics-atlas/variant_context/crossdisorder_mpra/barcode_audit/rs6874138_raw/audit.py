import csv,gzip,tarfile,collections,hashlib,json,math,re
from pathlib import Path
p=Path(__file__).resolve().parent;root=Path('/Volumes/Extreme SSD/brain6-crossdisorder-mpra');files={x:root/x for x in ('unique_bcdat.tsv.gz','GSE244011_RAW.tar','dnaMatrix.csv.gz','rnaMatrix.csv.gz')}
def sha(f):
 h=hashlib.sha256()
 with f.open('rb') as g:
  for b in iter(lambda:g.read(1048576),b''):h.update(b)
 return h.hexdigest()
prior=json.loads((p.parent/'coverage_receipt.json').read_text());actual={k:sha(v) for k,v in files.items()};print('hashes calculated',flush=True)
labels={}; mapped=collections.Counter()
with gzip.open(files['unique_bcdat.tsv.gz'],'rt',newline='') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  if r['rsid']!='rs6874138':continue
  v=r['variant'];m=re.match(r'^chr5_103899596_.*:T:C_([TC])$',v)
  if not m:continue
  allele='REF' if m[1]=='T' else 'ALT';code=r['barcode'];assert code not in labels or labels[code]==allele;labels[code]=allele;mapped[allele]+=1
assert mapped=={'REF':125,'ALT':131}
raw=[]
with tarfile.open(files['GSE244011_RAW.tar']) as t:
 members=t.getmembers();assert len(members)==44
 for i,item in enumerate(members):
  kind='DNA' if i<22 else 'RNA';assert kind in item.name.upper();counter=collections.Counter();unique={'REF':set(),'ALT':set()}
  with gzip.open(t.extractfile(item),'rt',newline='') as f:
   for r in csv.DictReader(f):
    allele=labels.get(r['barcode'])
    if allele:counter[allele]+=int(r['num_reads']);unique[allele].add(r['barcode'])
  raw.append({'kind':kind,'replicate_assumed':i%22+1,'source_file':item.name,'REF_raw':counter['REF'],'ALT_raw':counter['ALT'],'REF_barcodes':len(unique['REF']),'ALT_barcodes':len(unique['ALT'])})
print('raw files streamed',flush=True)
for kind in ('DNA','RNA'):
 with gzip.open(files[kind.lower()+'Matrix.csv.gz'],'rt',newline='') as f:
  matches=[r for r in csv.DictReader(f) if r['variant']=='chr5:103899596_T:C']
 assert len(matches)==1;row=matches[0]
 for rec in (x for x in raw if x['kind']==kind):
  for allele in ('REF','ALT'):
   col=f'{kind} Replicate {rec["replicate_assumed"]} {allele}';assert col in row;rec[allele+'_matrix']=int(row[col])
with (p/'counts.tsv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=raw[0],delimiter='\t');w.writeheader();w.writerows(raw)
def corr(x,y):
 mx=sum(x)/len(x);my=sum(y)/len(y);num=sum((a-mx)*(b-my) for a,b in zip(x,y));den=math.sqrt(sum((a-mx)**2 for a in x)*sum((b-my)**2 for b in y));return num/den if den else None
summaries=[]
for kind in ('DNA','RNA'):
 for allele in ('REF','ALT'):
  q=[r for r in raw if r['kind']==kind];x=[r[allele+'_raw'] for r in q];y=[r[allele+'_matrix'] for r in q];summaries.append({'kind':kind,'allele':allele,'raw_sum':sum(x),'matrix_sum':sum(y),'raw_to_matrix_sum_ratio':sum(x)/sum(y) if sum(y) else None,'assumed_order_pearson':corr(x,y),'equal_replicates':sum(a==b for a,b in zip(x,y)),'raw_greater_replicates':sum(a>b for a,b in zip(x,y))})
(p/'receipt.json').write_text(json.dumps({'source_sha256':actual,'mapped_labels':dict(mapped),'summary':summaries,'replicate_alignment':'assumed tar order; not independently verified'},indent=2)+'\n');print(json.dumps(summaries,indent=2))
