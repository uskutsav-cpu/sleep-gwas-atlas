#!/usr/bin/env python3
import csv,json,hashlib
import numpy as np,pandas as pd,pysam
url='https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr7.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz';v=pysam.VariantFile(url);panel=pd.read_csv('work/1kg_phase3_panel.tsv',sep='\t');samples=sorted(panel.loc[panel.super_pop.eq('EUR'),'sample']);pop=panel.set_index('sample').loc[samples,'pop'].to_numpy();wanted=[('rs1799252',113867592),('rs35624276',114040103),('rs12536335',114043159)];G=[];records=[]
for rs,pos in wanted:
 matches=[r for r in v.fetch('7',pos-1,pos) if r.pos==pos and rs in (r.id or '').split(';')];print(rs,'matching records',len(matches),flush=True)
 if len(matches)!=1:records.append({'rs':rs,'pos37':pos,'matching_records':len(matches)});continue
 r=matches[0];print('alleles',r.alleles,flush=True)
 if len(r.alleles)!=2:records.append({'rs':rs,'pos37':pos,'alleles':r.alleles,'status':'not biallelic'});continue
 vec=[]
 for s in samples:
  gt=r.samples[s].get('GT');vec.append(sum(g==1 for g in gt) if gt and len(gt)==2 and all(g is not None for g in gt) else np.nan)
 G.append(np.array(vec));records.append({'rs':rs,'pos37':pos,'alleles':r.alleles,'status':'biallelic','counted_allele':r.alleles[1]})
if len(G)!=3:print('not all variants matched, no pairwise correlations computed');open('outputs/stage6_device_sleep_foxp2_reference_ld.json','w').write(json.dumps({'variants':records,'source':url,'status':'incomplete'},indent=2));raise SystemExit(0)
rows=[]
for sub in ['EUR','FIN']:
 mask=np.ones(len(samples),dtype=bool) if sub=='EUR' else pop=='FIN'
 for k in [1,2]:
  ok=mask&np.isfinite(G[0])&np.isfinite(G[k]);r=float(np.corrcoef(G[0][ok],G[k][ok])[0,1]);rows.append({'reference':sub,'n':int(ok.sum()),'variant_1':wanted[0][0],'variant_2':wanted[k][0],'r_counted_alt':r,'r2':r*r,'allele_1':records[0]['counted_allele'],'allele_2':records[k]['counted_allele']})
pd.DataFrame(rows).to_csv('outputs/stage6_device_sleep_foxp2_reference_ld.csv',index=False);open('outputs/stage6_device_sleep_foxp2_reference_ld.json','w').write(json.dumps({'variants':records,'source':url,'status':'complete','note':'EUR includes FIN; signed r counts reference VCF alternate allele and not oriented to sleep-study effect.'},indent=2));print(pd.DataFrame(rows).to_string(index=False))
