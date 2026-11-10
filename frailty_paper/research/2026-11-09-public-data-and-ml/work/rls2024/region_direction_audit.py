import csv,json,math,hashlib,collections
from pathlib import Path
import numpy as np,pysam
# GRCh38 interval extends ~23 kb beyond each selected marker, deliberately set before retrieval.
chrom='2';start=66500000;end=66820000
urls={'akcimen':'https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90435001-GCST90436000/GCST90435387/harmonised/GCST90435387.h.tsv.gz','finngen':'https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_G6_RLS.gz'}
headers={'akcimen':['chrom','pos','ea','oa','beta','se','eaf','p','rsid','unused','hm_coordinate_conversion','hm_code','variant_id'],'finngen':['chrom','pos','ref','alt','rsid','nearest_genes','p','mlogp','beta','se','af_alt','af_alt_cases','af_alt_controls']}
store={};rawcounts={}
for source,url in urls.items():
 t=pysam.TabixFile(url); d=collections.defaultdict(list);rawcounts[source]=0
 for line in t.fetch(chrom,start-1,end):
  rawcounts[source]+=1;arr=line.rstrip('\r\n').split('\t');assert len(arr)==len(headers[source]);r=dict(zip(headers[source],arr))
  if source=='akcimen':
   if r['hm_code'] not in ('10','11','12','13','5','6'):continue
   try:ch,pos,ref,alt=r['variant_id'].split('_');pos=int(pos)
   except ValueError:continue
   if ch!=chrom or pos!=int(r['pos']) or {ref,alt}!={r['ea'],r['oa']}:continue
   beta=float(r['beta'])*(1 if r['ea']==alt else -1)
  else:ref=r['ref'];alt=r['alt'];pos=int(r['pos']);beta=float(r['beta'])
  if len(ref)!=1 or len(alt)!=1 or ref not in 'ACGT' or alt not in 'ACGT' or {ref,alt} in ({'A','T'},{'C','G'}):continue
  p=float(r['p']);se=float(r['se'])
  if not all(map(math.isfinite,[beta,p,se])) or p<0 or p>1 or se<=0:continue
  key=(pos,ref,alt);d[key].append((beta,se,p,r['rsid']))
 store[source]=d
out=[]
for key in sorted(set(store['akcimen'])&set(store['finngen'])):
 a=store['akcimen'][key];b=store['finngen'][key]
 if len(a)!=1 or len(b)!=1:continue
 x=a[0];y=b[0];out.append({'pos':key[0],'ref':key[1],'alt':key[2],'akcimen_rsid':x[3],'finngen_rsid':y[3],'akcimen_beta_alt':x[0],'akcimen_se':x[1],'akcimen_p':x[2],'finngen_beta_alt':y[0],'finngen_se':y[1],'finngen_p':y[2],'same_sign':int(x[0]*y[0]>0),'both_gws':int(x[2]<5e-8 and y[2]<5e-8)})
p=Path('outputs/stage6_rls_region_source_concordance.csv')
with p.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(out[0]));w.writeheader();w.writerows(out)
by=[]
for label,subset in [('all_matched',out),('both_p_below_5e-8',[r for r in out if r['both_gws']]),('within_25kb_rs113851554',[r for r in out if abs(r['pos']-66523432)<=25000 and r['both_gws']]),('within_25kb_rs78206187',[r for r in out if abs(r['pos']-66795102)<=25000 and r['both_gws']])]:
 by.append({'subset':label,'variants':len(subset),'same_sign':sum(r['same_sign'] for r in subset),'opposite_sign':sum(r['akcimen_beta_alt']*r['finngen_beta_alt']<0 for r in subset),'beta_correlation':float(np.corrcoef([r['akcimen_beta_alt'] for r in subset],[r['finngen_beta_alt'] for r in subset])[0,1]) if len(subset)>2 else None})
rec={'interval_GRCh38':f'{chrom}:{start}-{end}','sources':urls,'raw_rows':rawcounts,'filters':'unique exact chrom/pos/ref/alt non-palindromic SNVs; recognized harmonization codes; finite beta/se/p; excludes indels and duplicate exact records; does not rule out other alleles at a locus','subsets':by,'output_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'full_source_checksums_verified':False,'limitation':'Selected, LD-correlated region and two source p-values, not independent variants or a formal heterogeneity test; unknown study/relative overlap, imputation, phenotype and allele semantics.'};Path('outputs/stage6_rls_region_source_concordance_receipt.json').write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps(rec,indent=2))
