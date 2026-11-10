import csv
r=list(csv.DictReader(open('outputs/prototype_variant_five_trait_evidence.csv')))
traits=['FI','FFS','Jansen_insomnia','R13_insomnia','MVP_insomnia'];out=[]
for x in r:
 ps={k:float(x[k+'_p']) for k in traits}; betas={k:float(x[k+'_beta']) for k in traits}; signs={k:('positive' if v>0 else 'negative' if v<0 else 'zero') for k,v in betas.items()};
 out.append({'variant':x['variant'],'counted_allele':x['counted_allele'],'reference_multiallelic':x['reference_multiallelic'],'all_five_nominal_p_lt_0_05':all(p<.05 for p in ps.values()),'all_five_p_lt_0_05_div_15_descriptive':all(p<.05/15 for p in ps.values()),'all_five_signs_equal':len(set(signs.values()))==1,'largest_source_p':max(ps.values()),'largest_p_source':max(ps,key=ps.get),**{k+'_p':ps[k] for k in traits},**{k+'_beta_sign':signs[k] for k in traits},'status':'Post-selected descriptive screen; NOT multiplicity-adjusted confirmatory evidence, independent replication, or proof of common causal variant'})
with open('outputs/prototype_five_source_descriptive_triage.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=out[0]);w.writeheader();w.writerows(out)
assert len(out)==15 and sum(x['all_five_nominal_p_lt_0_05'] for x in out)==2 and sum(x['all_five_p_lt_0_05_div_15_descriptive'] for x in out)==0
print('Nominal five-source p<.05:',[(x['variant'],x['all_five_signs_equal']) for x in out if x['all_five_nominal_p_lt_0_05']]);print('Five source p<.05/15:',sum(x['all_five_p_lt_0_05_div_15_descriptive'] for x in out))
