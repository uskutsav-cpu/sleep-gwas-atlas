import csv,json,hashlib,collections,math,re
from pathlib import Path
import pysam
p=Path(__file__).resolve().parent;atlas=p.parents[1];base=Path('/Volumes/Brain6Work/xqtl-chr3');out=base/'ldz_readiness';out.mkdir(exist_ok=True)
old=atlas/'rs62264764_dense_qtl/multisignal_pilot';ref=list(csv.DictReader((old/'variants.tsv').open(),delimiter='\t'));assert len(ref)==1548
refkey=lambda r:(int(r['pos38']),r['ref'],r['alt'])
assert len({refkey(r) for r in ref})==len(ref)
def readrows(ctx):
 if ctx=='oli':
  x=pysam.TabixFile(str(base/'oli.bed.gz')); lines=x.fetch('chr3',117423157,118423158)
 else:lines=(base/'exc_region.tsv').open()
 d=collections.defaultdict(list)
 for line in lines:
  r=line.rstrip('\n').split('\t')
  if len(r)!=20 or r[9]!='ENSG00000239268':continue
  key=(int(r[2]),r[6],r[7]);m=re.search(r'(?:^|;)num_samples=(\d+)(?:;|$)',r[19]);
  if m:d[key].append(int(m.group(1)))
 if ctx=='oli':x.close()
 return d
counts={};fields=['variant','pos38','ref','alt','reference_alt_frequency','qtl_beta','qtl_se','qtl_alt_freq','num_samples','insomnia_beta','insomnia_se','adhd_beta','adhd_se']
for ctx in ('oli','exc'):
 samples=readrows(ctx);matched={}
 for t in ('insomnia','adhd'):
  f=base/('matched' if ctx=='oli' else 'exc_matched')/(ctx+'_'+t+'.tsv')
  rows=list(csv.DictReader(f.open(),delimiter='\t'));d=collections.defaultdict(list)
  for r in rows:d[(int(r['pos38']),r['ref'],r['alt'])].append(r)
  matched[t]=d
 rows=[];c=collections.Counter()
 for v in ref:
  key=refkey(v);a=matched['insomnia'].get(key,[]);b=matched['adhd'].get(key,[]);n=samples.get(key,[])
  if len(a)!=1 or len(b)!=1 or len(n)!=1:c['missing_or_duplicate']+=1;continue
  x,y=a[0],b[0];assert all(math.isclose(float(x[k]),float(y[k]),rel_tol=1e-12,abs_tol=1e-12) for k in ('qtl_beta','qtl_se','qtl_alt_freq'))
  for t,z in [('insomnia',x),('adhd',y)]:assert math.isclose(float(z['gwas_beta_alt']),float(v[t+'_beta']),rel_tol=1e-12,abs_tol=1e-12) and math.isclose(float(z['gwas_se']),float(v[t+'_se']),rel_tol=1e-12,abs_tol=1e-12)
  r={k:v[k] for k in ('variant','pos38','ref','alt','reference_alt_frequency','insomnia_beta','insomnia_se','adhd_beta','adhd_se')};r.update({'qtl_beta':x['qtl_beta'],'qtl_se':x['qtl_se'],'qtl_alt_freq':x['qtl_alt_freq'],'num_samples':n[0]});rows.append(r);c['matched']+=1
 assert len(rows)>1000
 dest=out/(ctx+'.tsv')
 with dest.open('w',newline='') as f:w=csv.DictWriter(f,fields,delimiter='\t');w.writeheader();w.writerows(rows)
 diffs=sorted(abs(float(r['qtl_alt_freq'])-float(r['reference_alt_frequency'])) for r in rows);ns=collections.Counter(r['num_samples'] for r in rows)
 counts[ctx]={'matched':len(rows),'counts':dict(c),'file':str(dest),'sha256':hashlib.sha256(dest.read_bytes()).hexdigest(),'num_samples_distribution':dict(ns),'frequency_abs_difference_median':diffs[len(diffs)//2],'frequency_abs_difference_max':max(diffs),'frequency_difference_above_0.1':sum(x>.1 for x in diffs),'frequency_difference_above_0.2':sum(x>.2 for x in diffs)}
 print(ctx,counts[ctx],flush=True)
(p/'matching_receipt.json').write_text(json.dumps({'cells':counts,'reference_marker_file_sha256':hashlib.sha256((old/'variants.tsv').read_bytes()).hexdigest(),'reference_ld_file_sha256':hashlib.sha256((old/'ld_matrix.tsv.gz').read_bytes()).hexdigest(),'limitations':'Matched small out-of-sample reference; matching frequency does not establish LD or source-model validity.'},indent=2)+'\n')
