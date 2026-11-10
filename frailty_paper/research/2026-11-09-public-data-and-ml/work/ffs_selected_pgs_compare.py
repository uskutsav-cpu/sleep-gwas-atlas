import csv,gzip,hashlib,json
from pathlib import Path
p=Path('outputs'); sf=Path('work/PGS005229_complete.txt.gz'); src=Path('/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/frailty_paper/data/gwas/fried_frailty_score/Frailty_GWAS_maf0.01_info0.9_hwe1e-6_FUMA.tsv')
def digest(f,alg):
 h=hashlib.new(alg)
 with f.open('rb') as q:
  for b in iter(lambda:q.read(4*1024*1024),b''):h.update(b)
 return h.hexdigest()
assert digest(sf,'md5')=='168e5c11469c0b27668a4c4d733276b1'
assert digest(src,'sha256')=='912e290b2064999c7147cbe3597f87cbcb14b25d387a8b6db2f25862e7274e55'
a={}
for label,file in [('15','unvalidated_15_variant_score_prototype.csv'),('22','unvalidated_22_shared_direction_score.csv')]:
 for r in csv.DictReader((p/file).open()):
  a.setdefault(r['variant'],[]).append((label,r))
ss={k:[] for k in a};ps={k:[] for k in a}
with src.open() as f:
 for r in csv.DictReader(f,delimiter='\t'):
  if r['SNP'] in a:ss[r['SNP']].append(r)
with gzip.open(sf,'rt') as f:
 for r in csv.DictReader((x for x in f if not x.startswith('#')),delimiter='\t'):
  if r['rsID'] in a:ps[r['rsID']].append(r)
out=[]
for k,items in a.items():
 for label,r in items:
  ch=r['chr'] if label=='15' else r['chr_GRCh37_assumed'];bp=r['position_GRCh37_inferred'] if label=='15' else r['position_GRCh37_assumed'];ca=r['counted_allele'];oa=r['other_allele']; match=[x for x in ss[k] if x['CHROM']==ch and x['POS']==bp and {x['A1'],x['A2']}=={ca,oa}];pp=[x for x in ps[k] if x['chr_name']==ch and x['chr_position']==bp];x=match[0] if len(match)==1 else {};y=pp[0] if len(pp)==1 else {}
  effect=y.get('effect_allele',''); other=y.get('other_allele','');out.append({'list':label,'variant':k,'chr':ch,'position':bp,'source_id_rows':len(ss[k]),'source_exact_pair_rows':len(match),'pgs_id_rows':len(ps[k]),'pgs_coord_rows':len(pp),'counted':ca,'other':oa,'source_A1':x.get('A1',''),'source_A2':x.get('A2',''),'source_beta':x.get('BETA',''),'prototype_ffs_beta':r.get('FFS_beta_weight_inferred',''),'pgs_effect':effect,'pgs_other':other,'pgs_weight':y.get('effect_weight',''),'pgs_effect_equals_source_A1':bool(x and y and effect==x['A1']),'pgs_other_valid_pair':bool(x and y and other in ('',x['A2']))})
with (p/'stage6_ffs_selected_pgs_compare.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=out[0].keys());w.writeheader();w.writerows(out)
print(json.dumps({'rows':len(out),'source_exact_pair':sum(x['source_exact_pair_rows']==1 for x in out),'pgs_coord':sum(x['pgs_coord_rows']==1 for x in out),'a1_matches':sum(x['pgs_effect_equals_source_A1'] for x in out),'other_valid':sum(x['pgs_other_valid_pair'] for x in out)},indent=2))
