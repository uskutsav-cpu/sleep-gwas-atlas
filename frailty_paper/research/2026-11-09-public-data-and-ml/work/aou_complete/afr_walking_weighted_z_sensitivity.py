import csv,json,math,hashlib
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import norm
O=Path('outputs');source=O/'stage6_aou_afr_walking_candidate_association.csv';df=pd.read_csv(source);assert len(df)==33 and (df.status=='matched').all() and (df.ancestry=='afr').all() and (df.sequencing_type=='genome').all()
sets=[('15_weighted',O/'unvalidated_15_variant_score_prototype.csv',O/'stage6_aou_weighted15_AFR661_ld.csv',[('FI','FI_beta_weight'),('FFS','FFS_beta_weight_inferred'),('insomnia','insomnia_beta_weight')]),('22_count',O/'unvalidated_22_shared_direction_score.csv',O/'stage6_aou_count22_AFR661_ld.csv',[('count','weight')])]
weighted22=pd.read_csv(O/'unvalidated_22_variant_score_prototype.csv').set_index('variant');res=[];components=[]
for name,weightsfile,ldfile,weights in sets:
 score=pd.read_csv(weightsfile);assert len(score)==(15 if name=='15_weighted' else 22)
 if name=='22_count':
  for label,col in [('FI','FI_beta_weight'),('FFS','FFS_beta_weight_inferred'),('insomnia','UKB_insomnia_beta_weight')]:
   vals=[]
   for r in score.itertuples(index=False):
    x=weighted22.loc[r.variant];assert {r.counted_allele,r.other_allele}=={x.counted_allele,x.other_allele}
    vals.append(float(x[col])*(1 if r.counted_allele==x.counted_allele else -1))
   score[label+'_oriented']=vals;weights.append((label,label+'_oriented'))
 ld=pd.read_csv(ldfile,index_col=0).loc[score.variant,score.variant].to_numpy();assert np.allclose(ld,ld.T) and np.linalg.eigvalsh(ld).min()>-1e-8
 for ph in ('NS_325.2',):
  sub=df[df.phenotype==ph].set_index('candidate').loc[score.variant];assert len(sub)==len(score)
  flip=[]
  for s,v in zip(score.itertuples(index=False),sub.itertuples()):
   assert {s.counted_allele,s.other_allele}=={v.ref,v.alt}
   flip.append(1 if s.counted_allele==v.alt else -1)
  beta=sub.beta.to_numpy(float);se=sub.se.to_numpy(float);p=sub.pvalue.to_numpy(float);freq=sub.af.to_numpy(float)
  assert ((p>0)&(p<=1)&(se>0)&(freq>0)&(freq<1)&(beta!=0)).all()
  z=norm.isf(p/2)*np.sign(beta)*np.array(flip);sd=np.sqrt(2*freq*(1-freq));precision=se*sd
  for label,col in weights:
   a=score[col].to_numpy(float)*sd;denom=np.sqrt(a@ld@a);assert denom>0
   Z=float(a@z/denom);drops=[]
   for i in range(len(a)):
    k=np.arange(len(a))!=i;drops.append(float(a[k]@z[k]/np.sqrt(a[k]@ld[np.ix_(k,k)]@a[k])))
   res.append({'score':name,'weights':label,'outcome':ph,'n_variants':len(a),'weighted_z_conditional':Z,'nominal_p_two_sided':float(2*norm.sf(abs(Z))),'leave_one_out_z_min':min(drops),'leave_one_out_z_max':max(drops),'se_times_genotype_sd_cv':float(precision.std(ddof=1)/precision.mean())})
   for i,s in enumerate(score.itertuples(index=False)):
    components.append({'score':name,'weights':label,'outcome':ph,'variant':s.variant,'conditional_counted_z':float(z[i]),'source_beta':float(beta[i]),'source_se':float(se[i]),'source_af':float(freq[i]),'hwe_genotype_sd':float(sd[i]),'weight_on_counted_allele':float(score[col].iloc[i]),'numerator_term':float(a[i]*z[i])})
p=O/'stage6_aou_afr_walking_weighted_z_sensitivity.csv';pd.DataFrame(res).to_csv(p,index=False);pd.DataFrame(components).to_csv(O/'stage6_aou_afr_walking_weighted_z_components.csv',index=False)
receipt={'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'output_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'reference':'1KG phase 3 AFR n=661; per-score counted allele orientation and within chromosome LD; cross-chromosome correlation zero','assumptions':'beta is ALT effect; biallelic diploid reversal; p to signed normal z; same samples/effective N; source AF is ALT frequency; HWE variance; AFR reference LD matches target Z covariance','limitations':'Postselected; case definitions differ; source overlap unverified; SPA and per-variant sample size/LD effects unaccounted for; reference samples may be related; FFS weights and multiallelic semantics uncertain; results are not individual-level prediction or validated score association.'}
(O/'stage6_aou_afr_walking_weighted_z_sensitivity_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(pd.DataFrame(res).to_string(index=False))
