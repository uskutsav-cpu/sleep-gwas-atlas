#!/usr/bin/env python3
"""Harmonize RBM6 model SNPs and six project GWAS regional caches."""
from pathlib import Path
import sqlite3,gzip,csv,json,math
from pyliftover import LiftOver
root=Path(__file__).resolve().parents[1];out=root/'outputs';work=root/'work';db=sqlite3.connect(f'file:{work}/PsychAD_snTWAS_EUR_Class_Mural.db?mode=ro',uri=True);weights=list(db.execute("select rsid,chr,weights.pos,ref_allele,eff_allele,weight from weights join extra using(gene) where weights.gene='ENSG00000004534'"));db.close(); assert len(weights)==31 and len(set(w[0] for w in weights))==31
wm={w[0]:w for w in weights}; lift=LiftOver(str(work/'hg19ToHg38.over.chain.gz'));specs={'Atkins_FI':('SNP','CHR','BP','A1','A2','BETA','SE',37),'Jansen_insomnia':('SNP','CHR','BP','A1','A2','BETA','SE',37),'Mak_FI':('SNP','CHR','POS','A1','A2','BETA','SE',37),'Mak_FP':('SNP','CHR','POS','A1','A2','BETA','SE',37),'Ye_FFS':('SNP','CHROM','POS','A1','A2','BETA','SE',37),'FinnGen_insomnia':('rsids','#chrom','pos','alt','ref','beta','sebeta',38),'FinnGen_HFRS_R12':('rs_id','chromosome','base_pair_location','effect_allele','other_allele','beta','standard_error',38)}
summary={}; audit=[]
for name,(idcol,chromcol,poscol,ecol,ocol,bcol,scol,build) in specs.items():
 found={k:[] for k in wm}
 with gzip.open(work/'stage3_region_cache'/f'{name}.regional.tsv.gz','rt') as f:
  for r in csv.DictReader(f,delimiter='\t'):
   if r[idcol] in found:found[r[idcol]].append(r)
 clean=[];counts={}
 for rs,w in wm.items():
  _,chrom,expected_pos,ref,eff,weight=w; records=found[rs];status='ok';beta=None;se=None
  if len(records)!=1:status='absent' if not records else 'duplicate_rsid'
  else:
   r=records[0]; c=r[chromcol].removeprefix('chr'); p=int(r[poscol]); e=r[ecol].upper();o=r[ocol].upper()
   mappings=lift.convert_coordinate('chr'+c,p-1) if build==37 else [('chr'+c,p-1,'+',0)]
   if len(mappings)!=1 or mappings[0][0]!='chr'+str(chrom) or mappings[0][1]+1!=expected_pos or mappings[0][2]!='+':status='coordinate_or_strand_mismatch'
   elif len(ref)!=1 or len(eff)!=1 or {ref,eff} not in ({'A','C'},{'A','G'},{'A','T'},{'C','G'},{'C','T'},{'G','T'}):status='invalid_model_alleles'
   elif {ref,eff} in ({'A','T'},{'C','G'}):status='palindromic_excluded'
   elif {e,o}!={ref,eff}:status='allele_mismatch'
   else:
    try:
     beta=float(r[bcol])*(1 if e==eff else -1);se=float(r[scol]);assert math.isfinite(beta) and math.isfinite(se) and se>0
    except:status='invalid_beta_or_se'
  counts[status]=counts.get(status,0)+1;audit.append({'trait':name,'rsid':rs,'model_pos_grch38':expected_pos,'model_other':ref,'model_effect':eff,'status':status,'gwas_beta_model_effect':beta if status=='ok' else '', 'gwas_se':se if status=='ok' else ''})
  if status=='ok':clean.append({'rsid':rs,'effect_allele':eff,'other_allele':ref,'beta':beta,'se':se})
 with open(work/f'stage6_rbm6_{name}_gwas.tsv','w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=['rsid','effect_allele','other_allele','beta','se'],delimiter='\t');w.writeheader();w.writerows(clean)
 summary[name]=counts
with open(out/'stage6_psychad_rbm6_gwas_harmonization.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(audit[0]));w.writeheader();w.writerows(audit)
(out/'stage6_psychad_rbm6_harmonization_receipt.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
