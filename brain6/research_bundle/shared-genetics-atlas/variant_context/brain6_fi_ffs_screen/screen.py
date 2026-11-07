import csv,gzip,hashlib,json,math,collections
from pathlib import Path
p=Path(__file__).resolve().parent;atlas=p.parents[1];pilot=atlas.parent/'brain6-alphagenome-pilot/results/jbrowse_expanded/expanded_variant_evidence.tsv'
fi=Path('/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/data/raw/frailty.txt.gz');ffs=Path('/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/frailty_paper/data/gwas/fried_frailty_score/Frailty_GWAS_maf0.01_info0.9_hwe1e-6_FUMA.tsv')
def sha(f):
 h=hashlib.sha256()
 with f.open('rb') as z:
  for b in iter(lambda:z.read(1048576),b''):h.update(b)
 return h.hexdigest()
assert sha(fi)=='a6702bc1c70bd1f2f3ec650d192dcdad72f7596bb716f8e6ef6c72359923e693'
assert sha(ffs)=='912e290b2064999c7147cbe3597f87cbcb14b25d387a8b6db2f25862e7274e55'
candidates=list(csv.DictReader(pilot.open(),delimiter='\t'));assert len(candidates)==425 and len({r['snp'] for r in candidates})==425
ids={r['snp'] for r in candidates};found={s:{k:[] for k in ids} for s in ('FI','FFS')};counts={}
for label,path in [('FI',fi),('FFS',ffs)]:
 n=0
 with (gzip.open(path,'rt') if label=='FI' else path.open()) as f:
  for r in csv.DictReader(f,delimiter='\t'):
   n+=1;key=r['variant_id'] if label=='FI' else r['SNP']
   if key in ids:found[label][key].append(r)
 counts[label]=n
out=[]
for c in candidates:
 chrom=c['chromosome'].removeprefix('chr');pos=c['pos_grch37_1based'];ea=c['source_a1'];oa=c['source_a2'];pal={ea,oa} in ({'A','T'},{'C','G'})
 for label in ('FI','FFS'):
  hits=found[label][c['snp']];valid=[]
  for r in hits:
   ch=r['chromosome'] if label=='FI' else r['CHROM'];bp=r['base_pair_location'] if label=='FI' else r['POS'];a=r['effect_allele'] if label=='FI' else r['A1'];b=r['other_allele'] if label=='FI' else r['A2']
   if ch.removeprefix('chr')==chrom and bp==pos and {a,b}=={ea,oa}:valid.append(r)
  rec={'variant':c['snp'],'region':c['candidate_locus_id'],'chr37':chrom,'pos37':pos,'effect_allele':ea,'other_allele':oa,'trait':label,'same_id_rows':len(hits),'exact_rows':len(valid),'status':'MISSING_ID' if not hits else 'NO_EXACT_MATCH' if not valid else 'DUPLICATE_EXACT' if len(valid)>1 else 'MATCHED'}
  if len(valid)==1:
   r=valid[0];a=r['effect_allele'] if label=='FI' else r['A1'];beta=r['beta'] if label=='FI' else r['BETA'];se=r['standard_error'] if label=='FI' else r['SE'];pv=r['p_value'] if label=='FI' else r['P']
   try:
    b=float(beta);s=float(se);v=float(pv);assert math.isfinite(b) and math.isfinite(s) and s>0 and math.isfinite(v) and 0<v<=1
    rec.update(p=v,se=s,source_effect_allele=a,beta_aligned='' if pal else (b if a==ea else -b),p_bonferroni_850=min(1,850*v),screen_hit=v<.05/850,orientation='PALINDROMIC_UNRESOLVED' if pal else 'SOURCE_A1_ASSUMED' if label=='FFS' else 'EFFECT_ALLELE')
   except (ValueError,AssertionError):rec['status']='INVALID_STATISTIC'
  out.append(rec)
keys=['variant','region','chr37','pos37','effect_allele','other_allele','trait','same_id_rows','exact_rows','status','p','se','source_effect_allele','beta_aligned','orientation','p_bonferroni_850','screen_hit']
with (p/'results.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,delimiter='\t',fieldnames=keys);w.writeheader();w.writerows(out)
summary={}
for t in ('FI','FFS'):
 rr=[x for x in out if x['trait']==t];hit=[x for x in rr if x.get('screen_hit')];summary[t]={'status_counts':dict(collections.Counter(x['status'] for x in rr)),'screen_hits':len(hit),'regions':dict(collections.Counter(x['region'] for x in hit)),'min_p':min((x['p'] for x in rr if 'p'in x),default=None)}
receipt={'source_rows':counts,'hashes':{'FI':sha(fi),'FFS':sha(ffs),'pilot':sha(pilot)},'candidate_count':425,'candidate_trait_rows':len(out),'family_slots':850,'threshold':.05/850,'summary':summary,'caveats':'postselected non-independent markers; no causal or corrected genome-wide discovery claim'}
(p/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
