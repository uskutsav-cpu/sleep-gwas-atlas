#!/usr/bin/env python3
import collections,csv,gzip,hashlib,json,math,statistics
from pathlib import Path
p=Path(__file__).resolve().parent;old=p.parents[2]/'rs62264764_dense_qtl/multisignal_pilot'
source=Path('/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/Sep1-local-dependencies/data/raw/adhd.txt.gz')
expected='c58a96031ec44b1edba81c91d603037a2efd03fcd946531424866e3f5f1d40c5'
h=hashlib.sha256()
with source.open('rb') as f:
 while b:=f.read(1024*1024):h.update(b)
assert h.hexdigest()==expected
def read(path):
 with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))
v={x['variant']:x for x in read(old/'variants.tsv')}; n={x['variant']:float(x['adhd_N_AUDIT']) for x in read(old/'sample_sizes.tsv')}
selected=[x for x in read(p/'repair_ldz_details.tsv') if x['trait']=='adhd'];variants={x['variant'] for x in selected}
assert len(variants)==1547
ids={v[x]['rsid'] for x in variants};assert len(ids)==1547
hits=collections.defaultdict(list)
with gzip.open(source,'rt') as f:
 head=f.readline().split(); ix={x:i for i,x in enumerate(head)}
 for line in f:
  row=line.split();rs=row[ix['SNP']]
  if rs in ids:hits[rs].append(dict(zip(head,row)))
assert len(hits)==1547 and all(len(x)==1 for x in hits.values())
data={};chars=collections.Counter();lengths=collections.Counter()
for variant in variants:
 a=v[variant];s=hits[a['rsid']][0]
 assert s['CHR']=='3' and s['BP']==a['pos37'] and {s['A1'],s['A2']}=={a['ref'],a['alt']}
 beta=math.log(float(s['OR']))*(1 if s['A1']==a['alt'] else -1)
 assert math.isclose(beta,float(a['adhd_beta']),abs_tol=1e-10,rel_tol=1e-10) and math.isclose(float(s['SE']),float(a['adhd_se']),abs_tol=1e-10,rel_tol=1e-10)
 case=float(s['Nca']);control=float(s['Nco']);eff=4*case*control/(case+control)
 assert math.isclose(eff,n[variant],abs_tol=1e-8)
 d=s['Direction'];chars.update(d);lengths[len(d)]+=1
 info=float(s['INFO']);assert math.isfinite(info)
 data[variant]={'mask':''.join('?' if c=='?' else '.' for c in d),'n':eff,'info':info}
assert set(chars)<=set('+-0?')
def summarize(vals):return {'count':len(vals),'median_neff':statistics.median(x['n'] for x in vals),'min_neff':min(x['n'] for x in vals),'max_neff':max(x['n'] for x in vals),'median_weighted_info':statistics.median(x['info'] for x in vals)}
out={}
for ref in ('eigenclip','ridge_0p001','1000g_eur'):
 rows=[r for r in selected if r['reference']==ref];assert len(rows)==1547
 groups=collections.defaultdict(list)
 for r in rows:
  x=data[r['variant']];flag='gt5' if abs(float(r['z_std_diff']))>5 else 'le5';groups[(x['mask'],flag)].append(x)
 out[ref]=[{'question_mark_mask':mask,'residual_group':flag,**summarize(vals)} for (mask,flag),vals in sorted(groups.items())]
receipt={'source_sha256':h.hexdigest(),'matched_variants':len(variants),'direction_character_counts':dict(chars),'direction_length_counts':dict(lengths),'reference_groups':out,'notes':'Question marks may mark missing cohort entries based on METAL source but release-specific encoding and cohort order not verified. Patterns replace all other characters with dots and do not retain effect signs. INFO is a cohort-weighted mean. Correlated variants, post-hoc grouping, and uncertain pairwise sample overlap prevent causal interpretation.'}
(p/'adhd_direction_audit.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
