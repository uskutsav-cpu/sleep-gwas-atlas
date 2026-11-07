#!/usr/bin/env python3
import csv,collections,json,math,hashlib
from pathlib import Path
root=Path('outputs/shared-genetics-atlas/marginal_priority');out=root/'hfrs_followup'
source=Path('/Volumes/Extreme SSD/sleep-frailty-full-autonomous-v5-20261009/frailty_paper/full_autonomous_research_v5/genetics/new_public_validation/filtered/HFRS_HM3.tsv')
with (root/'marginal_threshold_crossings.tsv').open() as f:sel=list(csv.DictReader(f,delimiter='\t'))
assert len(sel)==368 and len({r['SNP'] for r in sel})==len(sel)
wanted={r['SNP'] for r in sel};matched=collections.defaultdict(list);h=hashlib.sha256()
with source.open('rb') as f:
 while b:=f.read(4*1024*1024):h.update(b)
with source.open() as f:
 for r in csv.DictReader(f,delimiter='\t'):
  if r['SNP'] in wanted:matched[r['SNP']].append(r)
cut=.05/(4*1138262);status=collections.Counter();patterns=collections.Counter();results=[]
for a in sel:
 r=dict(a);m=matched[a['SNP']];r.update(hfrs_status='',hfrs_A1='',hfrs_A2='',hfrs_beta_A1='',hfrs_beta_aligned='',hfrs_p='',four_trait_pattern='')
 if len(m)!=1:r['hfrs_status']='MISSING' if not m else 'DUPLICATE'
 else:
  b=m[0];aa=(a['A1'],a['A2']);bb=(b['A1'],b['A2'])
  if a['chr37']!=b['CHR'] or len(set(aa))!=2 or any(len(x)!=1 or x not in 'ACGT' for x in aa+bb) or set(aa)!=set(bb):r['hfrs_status']='CHR_OR_ALLELE_MISMATCH'
  else:
   try:
    p=float(b['P']); beta=float(b['BETA']);assert math.isfinite(beta) and math.isfinite(p) and 0<p<=1
    r.update(hfrs_status='MATCHED',hfrs_A1=bb[0],hfrs_A2=bb[1],hfrs_beta_A1=beta,hfrs_beta_aligned=beta if aa==bb else -beta,hfrs_p=p)
    names=['insomnia','FI','FFS','HFRS'];pvals=[float(a['p_insomnia_upper']),float(a['p_FI_upper']),float(a['p_FFS_upper']),p]
    pattern='+'.join(name for name,value in zip(names,pvals) if value<=cut) or 'NONE';r['four_trait_pattern']=pattern;patterns[pattern]+=1
   except (ValueError,AssertionError):r['hfrs_status']='INVALID'
 status[r['hfrs_status']]+=1;results.append(r)
with (out/'selected_with_hfrs.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(results[0]),delimiter='\t');w.writeheader();w.writerows(results)
summary={'selection_n':len(sel),'hfrs_source':str(source),'hfrs_source_sha256':h.hexdigest(),'selected_source_sha256':hashlib.sha256((root/'marginal_threshold_crossings.tsv').read_bytes()).hexdigest(),'threshold_four_trait':cut,'status':dict(status),'patterns_on_matched_selected':dict(patterns),'note':'Selected three-trait SNPs only; does not survey HFRS-only variants, test trait specificity or count independent loci.'}
(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
