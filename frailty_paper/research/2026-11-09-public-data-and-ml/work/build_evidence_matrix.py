import csv,math,collections
from pathlib import Path
proto=list(csv.DictReader(open('outputs/unvalidated_15_variant_score_prototype.csv')))
main={x['candidate']:x for x in csv.DictReader(open('outputs/published_candidate_marginal_lookup.csv'))}
external={(x['variant'],x['source']):x for x in csv.DictReader(open('outputs/prototype_15_variants_full_external_insomnia_lookup.csv'))}
rows=[]
for x in proto:
 variant=x['variant']; m=main[variant]; o={'variant':variant,'chromosome_GRCh37':x['chr'],'position_GRCh37':x['position_GRCh37_inferred'],'counted_allele':x['counted_allele'],'other_allele':x['other_allele'],'reference_multiallelic':x['GRCh37_reference_multiallelic'],'published_cross_trait_rank':x['song_cross_trait_rank']}
 for label,beta,se,p in [('FI',x['FI_beta_weight'],m['FI_SE'],m['FI_P']),('FFS',x['FFS_beta_weight_inferred'],m['FFS_SE'],m['FFS_P']),('Jansen_insomnia',x['insomnia_beta_weight'],m['insomnia_SE'],m['insomnia_P'])]:
  beta=float(beta);se=float(se);o.update({label+'_beta':beta,label+'_se':se,label+'_ci95_low':beta-1.959963984540054*se,label+'_ci95_high':beta+1.959963984540054*se,label+'_p':p})
 for label in ('R13','MVP'):
  z=external[(variant,label)];beta=float(z['aligned_beta_to_candidate_allele']); se=None if z['source_SE'] in ('#NA','','NA') else float(z['source_SE']);o.update({label+'_insomnia_beta':beta,label+'_insomnia_se':se if se is not None else '',label+'_insomnia_ci95_low':beta-1.959963984540054*se if se is not None else '',label+'_insomnia_ci95_high':beta+1.959963984540054*se if se is not None else '',label+'_insomnia_p':z['source_p']})
 o['FI_and_FFS_nominal_p_below_0.05']=str(float(o['FI_p'])<.05 and float(o['FFS_p'])<.05)
 o['FI_and_FFS_aligned_signs_same']=str(float(o['FI_beta'])*float(o['FFS_beta'])>0)
 o['status']='Exploratory post-selected candidates; CIs unadjusted; FI/FFS sample overlap unknown; FFS allele convention inferred; multiallelic contrasts unresolved; no risk probability or evidence of shared causality'
 rows.append(o)
assert len(rows)==15
with open('outputs/prototype_variant_five_trait_evidence.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
print('rows',len(rows),'FI/FFS nominal p<.05',sum(x['FI_and_FFS_nominal_p_below_0.05']=='True' for x in rows),'same sign among nominal',sum(x['FI_and_FFS_nominal_p_below_0.05']=='True' and x['FI_and_FFS_aligned_signs_same']=='True' for x in rows),'MVP CIs',sum(x['MVP_insomnia_se']!='' for x in rows))
for x in rows:
 if x['FI_and_FFS_nominal_p_below_0.05']=='True':print(x['variant'],x['FI_p'],x['FFS_p'],x['Jansen_insomnia_p'],x['R13_insomnia_p'],x['MVP_insomnia_p'])
