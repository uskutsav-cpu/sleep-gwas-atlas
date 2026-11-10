#!/usr/bin/env python3
"""Exploratory score-weighted marginal Z test; not individual-level score validation."""
import csv,json,hashlib
from pathlib import Path
import numpy as np,pandas as pd,pysam
from scipy.stats import norm
O=Path('outputs');W=Path('work');score=pd.read_csv(O/'unvalidated_15_variant_score_prototype.csv');g=pd.read_csv(O/'candidate_finngen_hfrs_full_gwas_lookup.csv');d=score.merge(g,left_on='variant',right_on='variant',validate='one_to_one');assert len(d)==15
assert (d.counted_allele==d.candidate_counted_allele).all() and (d.other_allele==d.candidate_other_allele).all();assert d.hfrs_beta_on_counted_allele.notna().all() and (d.hfrs_se>0).all()
panel=pd.read_csv(W/'1kg_phase3_panel.tsv',sep='\t');pops=dict(zip(panel['sample'],panel['pop']));samples=sorted(panel.loc[panel['super_pop'].eq('EUR'),'sample']);assert len(samples)==503
url='https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz';files={};G=[];checks=[]
for row in d.itertuples():
 chrom=str(row.chr);pos=int(row.position_GRCh37_inferred);a=row.counted_allele;b=row.other_allele
 if chrom not in files:files[chrom]=pysam.VariantFile(url.format(chrom=chrom))
 f=files[chrom];matches=[x for x in f.fetch(chrom,pos-1,pos) if x.pos==pos and x.id==row.variant];assert len(matches)==1
 record=matches[0];assert len(record.alleles)==2 and set(record.alleles)=={a,b};i=record.alleles.index(a);vec=[]
 for sample in samples:
  gt=record.samples[sample].get('GT');assert gt is not None and len(gt)==2 and all(v is not None for v in gt);vec.append(sum(v==i for v in gt))
 G.append(vec);checks.append({'variant':row.variant,'alleles':list(record.alleles),'EUR_called':len(vec)})
 print(row.variant,'matched',flush=True)
G=np.asarray(G,dtype=float).T;np.savez_compressed(W/'stage6_hfrs_sparse_score_EUR503_genotypes.npz',genotypes=G,variants=d.variant.to_numpy(),samples=np.array(samples));R=np.corrcoef(G,rowvar=False);ch=d['chr'].astype(str).to_numpy();R=np.where(ch[:,None]==ch[None,:],R,0.0);assert np.linalg.eigvalsh(R).min()>-1e-8
# Target GWAS frequency is for its source effect allele, mapped to counted allele.
f=np.where(d.hfrs_source_effect_allele==d.counted_allele,d.hfrs_source_effect_allele_frequency,1-d.hfrs_source_effect_allele_frequency);assert np.all((f>0)&(f<1));std=np.sqrt(2*f*(1-f));z=(d.hfrs_beta_on_counted_allele/d.hfrs_se).to_numpy();
rows=[];detail=[]
for weight in ['FI_beta_weight','FFS_beta_weight_inferred','insomnia_beta_weight']:
 w=d[weight].to_numpy()*std
 for omit in [None,*d.variant]:
  keep=np.ones(len(d),dtype=bool) if omit is None else (d.variant!=omit).to_numpy();a=w[keep];zz=z[keep];RR=R[np.ix_(keep,keep)];denom=float(np.sqrt(a@RR@a));stat=float(a@zz/denom);rows.append({'weight':weight,'omitted_variant':omit or 'none','n':sum(keep),'weighted_z':stat,'two_sided_p':2*norm.sf(abs(stat)),'sd_assuming_independent_markers':float(np.sqrt(a@a)),'sd_with_reference_ld':denom})
 for i,variant in enumerate(d.variant):detail.append({'weight':weight,'variant':variant,'target_z':z[i],'target_counted_allele_freq':f[i],'genotype_sd_hwe':std[i],'weight_per_counted_allele':d[weight][i],'weighted_z_numerator_component':w[i]*z[i]})
pd.DataFrame(rows).to_csv(O/'stage6_hfrs_sparse_score_summary_test_EUR503.csv',index=False);pd.DataFrame(detail).to_csv(O/'stage6_hfrs_sparse_score_summary_components_EUR503.csv',index=False);pd.DataFrame(R,index=d.variant,columns=d.variant).to_csv(O/'stage6_hfrs_sparse_score_EUR503_ld.csv');receipt={'source_reference':url,'samples':503,'reference_population':'1000 Genomes phase 3 EUR (five subpopulations)','markers':checks,'cross_chromosome_ld_assumed_zero':True,'within_chromosome_ld':'sample dosage correlations','genotype_sd':'Hardy-Weinberg approximation sqrt(2 p (1-p)) using HFRS reported counted allele frequency','gwass':'GCST90707291 exact lookup','weight_sha256':hashlib.sha256((O/'unvalidated_15_variant_score_prototype.csv').read_bytes()).hexdigest(),'limitation':'null Z-correlation approximated by reference LD, assumes compatible study samples and errors, approximate common effective sample and additive model, not held-out individual-level score prediction'};(O/'stage6_hfrs_sparse_score_summary_receipt_EUR503.json').write_text(json.dumps(receipt,indent=2)+'\n');print(pd.DataFrame(rows).query("omitted_variant == 'none'").to_string(index=False))
