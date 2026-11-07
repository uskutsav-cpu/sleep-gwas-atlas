#!/usr/bin/env python3
import collections,csv,gzip,hashlib,json,math,statistics
from pathlib import Path
p=Path(__file__).resolve().parent; old=p.parents[2]/'rs62264764_dense_qtl/multisignal_pilot'
source=Path('/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/Sep1-local-dependencies/data/raw/adhd.txt.gz')
expected='c58a96031ec44b1edba81c91d603037a2efd03fcd946531424866e3f5f1d40c5'
h=hashlib.sha256()
with source.open('rb') as f:
 while b:=f.read(1024*1024):h.update(b)
assert h.hexdigest()==expected
def tab(path):
 with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))
v={x['variant']:x for x in tab(old/'variants.tsv')}
n={x['variant']:float(x['adhd_N_AUDIT']) for x in tab(old/'sample_sizes.tsv')}
g={x['variant']:x for x in tab(p/'matched_snps.tsv')}
diag=[x for x in tab(p/'repair_ldz_details.tsv') if x['reference']=='eigenclip' and x['trait']=='adhd']
top=sorted(diag,key=lambda x:abs(float(x['z_std_diff'])),reverse=True)[:10]
assert len(top)==10 and len({x['variant'] for x in top})==10
ids={v[x['variant']]['rsid'] for x in top}; assert len(ids)==10
hits=collections.defaultdict(list)
with gzip.open(source,'rt') as f:
 head=f.readline().split(); ix={x:i for i,x in enumerate(head)}
 for line in f:
  row=line.split(); rs=row[ix['SNP']]
  if rs in ids:hits[rs].append(dict(zip(head,row)))
assert set(hits)==ids and all(len(rows)==1 for rows in hits.values())
ds_old=[];ds_new=[];neffs=[];out=[]
for x in top:
 variant=x['variant']; a=v[variant]; b=g[variant]; s=hits[a['rsid']][0]
 assert s['CHR']=='3' and s['BP']==a['pos37'] and {s['A1'],s['A2']}=={a['ref'],a['alt']}
 beta=math.log(float(s['OR']))*(1 if s['A1']==a['alt'] else -1)
 assert math.isclose(beta,float(a['adhd_beta']),rel_tol=1e-10,abs_tol=1e-10)
 assert math.isclose(float(s['SE']),float(a['adhd_se']),rel_tol=1e-10,abs_tol=1e-10)
 f=float(s['FRQ_U_186843']); source_alt=f if s['A1']==a['alt'] else 1-f
 gctb_alt=float(b['gctb_a1freq']) if b['gctb_a1']==a['alt'] else 1-float(b['gctb_a1freq'])
 ds_old.append(abs(source_alt-float(a['reference_alt_frequency'])));ds_new.append(abs(source_alt-gctb_alt))
 case=float(s['Nca']);control=float(s['Nco']);eff=4*case*control/(case+control)
 assert math.isclose(eff,n[variant],abs_tol=1e-8);neffs.append(eff)
 out.append({'variant':variant,'rsid':a['rsid'],'std_residual':float(x['z_std_diff'])})
def describe(x):return {'median':statistics.median(x),'maximum':max(x),'above_0_05':sum(y>.05 for y in x),'above_0_1':sum(y>.1 for y in x)}
result={'source_sha256':h.hexdigest(),'selection':'top 10 absolute ADHD residuals in GCTB eigenclip with N 128214','unique_matches_and_aligned_beta_se':len(top),'top_variants':out,'source_control_vs_1000g_alt_frequency_difference':describe(ds_old),'source_control_vs_gctb_alt_frequency_difference':describe(ds_new),'audited_variant_neff_min':min(neffs),'audited_variant_neff_max':max(neffs),'limitations':'Post-selected sample; control frequencies can differ from reference populations, and cannot validate LD or causal interpretations. No original source rows redistributed.'}
(p/'adhd_source_audit.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
