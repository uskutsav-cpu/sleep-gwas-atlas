import csv,pathlib,collections,math,json,hashlib
root=pathlib.Path('/Volumes/Extreme SSD/Codex Archive/mak2023_same_cohort');rows=list(csv.DictReader(open('outputs/expanded_22_candidate_risk_gene_table.csv')));wanted={r['variant']:r for r in rows};found=collections.defaultdict(list);n=0
with (root/'fi.fastGWA').open() as fi,(root/'fp.fastGWA').open() as fp:
 a=csv.DictReader(fi,delimiter='\t');b=csv.DictReader(fp,delimiter='\t')
 for u,v in zip(a,b):
  n+=1;assert (u['SNP'],u['CHR'],u['POS'])==(v['SNP'],v['CHR'],v['POS']);rs=u['SNP']
  if rs in wanted:found[rs].append((u,v))
result=[]
for x in rows:
 rs=x['variant'];same=[(u,v) for u,v in found[rs] if u['CHR']==x['GRCh37_chr_assumed'] and u['POS']==x['GRCh37_pos_assumed'] and {u['A1'],u['A2']}=={x['counted_allele'],x['other_allele']} and {v['A1'],v['A2']}=={x['counted_allele'],x['other_allele']}]
 d={'variant':rs,'counted_allele':x['counted_allele'],'other_allele':x['other_allele'],'positional_gene':x['nearest_gene_body_positional'],'source_id_hits':len(found[rs]),'exact_hits':len(same),'prior_FI_beta':x['FI_beta_aligned'],'prior_FI_p':x['FI_p'],'prior_FFS_beta':x['FFS_beta_aligned'],'prior_FFS_p':x['FFS_p']}
 if len(same)==1:
  u,v=same[0]
  for tag,r in (('mak_FI',u),('mak_FP',v)):
   sign=1 if r['A1']==x['counted_allele'] else -1;beta=sign*float(r['BETA']);se=float(r['SE']);assert math.isfinite(beta) and math.isfinite(se) and se>0
   d.update({tag+'_beta':beta,tag+'_se':se,tag+'_ci_low':beta-1.95996398454*se,tag+'_ci_high':beta+1.95996398454*se,tag+'_p':float(r['P']),tag+'_N':r['N'],tag+'_A1':r['A1'],tag+'_AF1':r['AF1']})
  d.update(FI_same_sign= d['mak_FI_beta']*float(x['FI_beta_aligned'])>0,FP_same_sign=d['mak_FP_beta']*float(x['FFS_beta_aligned'])>0,mak_both_nominal=d['mak_FI_p']<.05 and d['mak_FP_p']<.05)
 result.append(d)
out=pathlib.Path('outputs/mak2023_22_candidate_lookup.csv');fields=list(dict.fromkeys(k for d in result for k in d))
with out.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(result)
valid=[x for x in result if x['exact_hits']==1];rec={'source_paired_rows':n,'unique_matches':len(valid),'FI_same_sign':sum(x['FI_same_sign'] for x in valid),'FP_same_sign':sum(x['FP_same_sign'] for x in valid),'both_nominal':sum(x['mak_both_nominal'] for x in valid),'both_nominal_same_sign':[x['variant'] for x in valid if x['mak_both_nominal'] and x['FI_same_sign'] and x['FP_same_sign']],'output_sha256':hashlib.sha256(out.read_bytes()).hexdigest(),'source_MD5_verified_in_earlier_analysis':{'fi':'5cd10928bde27e48b82d38f097696473','fp':'3ec3de20038f3022eaa9163347b0e194'}}
pathlib.Path('outputs/mak2023_22_candidate_lookup.json').write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps(rec,indent=2))
