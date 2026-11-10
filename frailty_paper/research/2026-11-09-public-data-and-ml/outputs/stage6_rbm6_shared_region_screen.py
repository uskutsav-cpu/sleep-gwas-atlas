#!/usr/bin/env python3
"""Descriptive exact-allele regional comparison; this is not colocalization."""
from pathlib import Path
import gzip,csv,json,collections,math
import numpy as np,pandas as pd
from pyliftover import LiftOver
root=Path(__file__).resolve().parents[1];work=root/'work';out=root/'outputs';lo,hi=49008613,51008613;lift=LiftOver(str(work/'hg19ToHg38.over.chain.gz'))
specs={'FinnGen_insomnia':('stage3_region_cache/FinnGen_insomnia.regional.tsv.gz','rsids','#chrom','pos','alt','ref','beta','sebeta','pval',38),'FinnGen_HFRS':('stage6_rbm6_hfrs_2mb.tsv.gz','rs_id','chromosome','base_pair_location','effect_allele','other_allele','beta','standard_error','p_value',38),'Jansen_insomnia':('stage3_region_cache/Jansen_insomnia.regional.tsv.gz','SNP','CHR','BP','A1','A2','BETA','SE','P',37),'Mak_FI':('stage3_region_cache/Mak_FI.regional.tsv.gz','SNP','CHR','POS','A1','A2','BETA','SE','P',37),'Ye_FFS':('stage3_region_cache/Ye_FFS.regional.tsv.gz','SNP','CHROM','POS','A1','A2','BETA','SE','P',37)}
allrows={};audit={}
for name,(file,idc,cc,pc,ac,bc,betac,sec,pcol,build) in specs.items():
 rows=[];count=collections.Counter()
 with gzip.open(work/file,'rt') as f:
  for r in csv.DictReader(f,delimiter='\t'):
   if r[cc].removeprefix('chr')!='3':continue
   try:p=int(r[pc]);a=r[ac].upper();b=r[bc].upper()
   except:continue
   if build==37:
    # Restrict rough pre-liftover span to avoid mapping all chromosome 3 cached variants.
    if not 48500000<=p<=51500000:continue
    m=lift.convert_coordinate('chr3',p-1)
    if len(m)!=1 or m[0][0]!='chr3' or m[0][2]!='+':count['bad_liftover']+=1;continue
    p=m[0][1]+1
   if not lo<=p<=hi:continue
   if len(a)!=1 or len(b)!=1 or {a,b} not in ({'A','C'},{'A','G'},{'A','T'},{'C','G'},{'C','T'},{'G','T'}):count['not_biallelic_snp']+=1;continue
   if {a,b} in ({'A','T'},{'C','G'}):count['palindromic']+=1;continue
   try:beta=float(r[betac]);se=float(r[sec]);pval=float(r[pcol]);assert all(map(math.isfinite,(beta,se,pval))) and se>0 and 0<=pval<=1
   except:count['invalid_statistics']+=1;continue
   rs=r[idc]
   if not rs.startswith('rs') or not rs[2:].isdigit():count['invalid_or_multi_rsid']+=1;continue
   rows.append({'pos':p,'rsid':rs,'effect':a,'other':b,'pair':''.join(sorted([a,b])),'beta':beta,'se':se,'p':pval})
 d=pd.DataFrame(rows);count['valid_before_duplicates']=len(d);dups=d.duplicated(['pos','rsid','pair'],keep=False);count['duplicate_rows_removed']=int(dups.sum());d=d[~dups].copy();audit[name]=dict(count);allrows[name]=d
# FinnGen R13 ref allele is source reference; use it for a common GRCh38 ALT orientation.
ref=allrows['FinnGen_insomnia'].copy();ref=ref.rename(columns={'effect':'alt','other':'ref','beta':'FinnGen_insomnia_beta_alt','se':'FinnGen_insomnia_se','p':'FinnGen_insomnia_p'}); cols=['pos','rsid','pair','ref','alt','FinnGen_insomnia_beta_alt','FinnGen_insomnia_se','FinnGen_insomnia_p'];merged=ref[cols]
for name in ['FinnGen_HFRS','Jansen_insomnia','Mak_FI','Ye_FFS']:
 d=allrows[name]; merged=merged.merge(d,on=['pos','rsid','pair'],how='inner',validate='one_to_one');merged[name+'_beta_alt']=np.where(merged.effect==merged.alt,merged.beta,-merged.beta);merged[name+'_se']=merged.se;merged[name+'_p']=merged.p;merged=merged.drop(columns=['effect','other','beta','se','p'])
merged=merged.sort_values('pos');merged.to_csv(out/'stage6_rbm6_five_gwas_common_region.csv.gz',index=False,compression='gzip')
results={'region':'GRCh38 chr3:49008613-51008613','filters':audit,'all_five_shared_snps':len(merged),'comparisons':{}}
for sleep,frail in [('Jansen_insomnia','FinnGen_HFRS'),('Jansen_insomnia','Mak_FI'),('FinnGen_insomnia','FinnGen_HFRS'),('FinnGen_insomnia','Ye_FFS')]:
 df=merged;zs=df[sleep+'_beta_alt']/df[sleep+'_se'];zf=df[frail+'_beta_alt']/df[frail+'_se'];i=df[sleep+'_p'].idxmin();j=df[frail+'_p'].idxmin()
 def item(k):
  r=df.loc[k];return {'rsid':r.rsid,'pos':int(r.pos),'ref':r.ref,'alt':r.alt,'sleep_p':float(r[sleep+'_p']),'frail_p':float(r[frail+'_p'])}
 results['comparisons'][sleep+'_vs_'+frail]={'n':len(df),'signed_z_pearson_descriptive':float(zs.corr(zf)),'sleep_best_shared_variant':item(i),'frailty_best_shared_variant':item(j)}
(out/'stage6_rbm6_five_gwas_region_summary.json').write_text(json.dumps(results,indent=2)+'\n');print(json.dumps(results,indent=2))
