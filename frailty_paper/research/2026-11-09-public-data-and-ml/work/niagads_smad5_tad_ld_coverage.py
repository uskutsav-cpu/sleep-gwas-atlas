#!/usr/bin/env python3
"""Audit available EUR reference coverage before fitting a wider-window model."""
import csv, hashlib, json, math
from pathlib import Path
import numpy as np

root=Path(__file__).resolve().parent.parent;out=root/'outputs'
ref=root/'work/stage4_magma_reference'; n=503
traits={}
for name in ('Mak_FI','Mak_FP','FinnGen_insomnia'):
 with (out/f'stage6_niagads_smad5_exc_{name}_aligned_tad.csv').open() as f:
  traits[name]=[r for r in csv.DictReader(f) if abs(float(r['qtl_alt_af'])-float(r['gwas_af_alt']))<=.1]
all_rows={}
for rows in traits.values():
 for r in rows:
  _,pos,refallele,alt=r['variant_GRCh38'].split(':')
  key=(int(r['gwas_pos37']),frozenset([refallele,alt]))
  all_rows.setdefault(key,set()).add(r['variant_GRCh38'])
matches={}
with (ref/'g1000_eur.bim').open() as f:
 for i,line in enumerate(f):
  z=line.split()
  if z[0]!='5': continue
  key=(int(z[3]),frozenset((z[4],z[5])))
  if key in all_rows: matches.setdefault(key,[]).append((i,z[4],z[5]))
fam=(ref/'g1000_eur.fam').read_text().splitlines();assert len(fam)==n
bed=np.memmap(ref/'g1000_eur.bed',mode='r',dtype='uint8');assert bytes(bed[:3])==bytes([108,27,1]); step=(n+3)//4
retained={}
for key,entries in matches.items():
 if len(entries)!=1 or len(all_rows[key])!=1:continue
 var=next(iter(all_rows[key]));alt=var.split(':')[-1];idx,a1,a2=entries[0]
 raw=bed[3+idx*step:3+(idx+1)*step]
 codes=((raw[:,None]>>np.array([0,2,4,6],dtype=np.uint8))&3).reshape(-1)[:n]
 g=np.array([2.,np.nan,1.,0.])[codes]; g=g if a1==alt else 2-g
 af=float(np.nanmean(g)/2)
 if not np.isfinite(af) or af<=0 or af>=1:continue
 retained[var]={'index':idx,'reference_alt_frequency':af,'missing':int(np.isnan(g).sum())}
def bf(beta,se,sd):
 w=sd*sd;v=se*se;rat=w/(w+v);return .5*(math.log1p(-rat)+rat*(beta/se)**2)
def support(weights,mask):
 top=max(weights);a=np.exp(np.array(weights)-top);return float(a[mask].sum()/a.sum())
res={}
for name,rows in traits.items():
 mask=np.array([r['variant_GRCh38'] in retained and abs(retained[r['variant_GRCh38']]['reference_alt_frequency']-float(r['qtl_alt_af']))<=.1 for r in rows])
 q=[bf(float(r['qtl_beta_alt']),float(r['qtl_se']),.15) for r in rows]
 sd={'Mak_FI':1.1004195728364333,'Mak_FP':.1232461359076205,'FinnGen_insomnia':.2}[name]
 g=[bf(float(r['gwas_beta_alt']),float(r['gwas_se']),sd) for r in rows]
 top=rows[int(np.argmax(g))]['variant_GRCh38']
 res[name]={'input_n':len(rows),'matched_reference_n':int(mask.sum()),
            'qtl_abf_support_fraction':support(q,mask),'gwas_abf_support_fraction':support(g,mask),
            'gwas_top_variant':top,'gwas_top_retained':bool(mask[int(np.argmax(g))]),
            'expression_lead_retained':any(r['variant_GRCh38']=='chr5:136206558:C:T' and m for r,m in zip(rows,mask)),
            'candidate_retained':any(r['variant_GRCh38']=='chr5:136202653:A:C' and m for r,m in zip(rows,mask))}
 if name=='FinnGen_insomnia':
  with (root/'work/niagads_smad5_tad_insomnia_ld_inputs.csv').open('w',newline='') as f:
   writer=csv.DictWriter(f,fieldnames=list(rows[0])+['reference_bim_index','reference_alt_frequency']);writer.writeheader()
   for row,keep in zip(rows,mask):
    if keep:
     item=retained[row['variant_GRCh38']];writer.writerow({**row,'reference_bim_index':item['index'],'reference_alt_frequency':item['reference_alt_frequency']})
receipt={'reference_n':n,'reference':'1000 Genomes EUR local PLINK reference; 503 samples',
         'reference_bim_sha256':hashlib.sha256((ref/'g1000_eur.bim').read_bytes()).hexdigest(),
         'results':res,'limitations':'Support fractions conditional on matched filtered TAD set and assumed priors, not complete GWAS coverage. Reference is not target LD; no score validation.'}
(out/'stage6_niagads_smad5_tad_ld_coverage.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(res,indent=2))
