import csv,math,hashlib,pathlib,json
p=pathlib.Path('outputs'); vs=['rs12536335','rs1909650']
def read(name):return {r['variant']:r for r in csv.DictReader((p/name).open())}
five=read('prototype_variant_five_trait_evidence.csv'); h=read('candidate_finngen_hfrs_full_gwas_lookup.csv'); tg=read('candidate_twingene_fi_cohort_lookup.csv'); raw=[r for r in csv.DictReader(open('work/raw_external_prototype_matches.csv')) if r['source']=='MVP']; raw={r['variant']:r for r in raw};out=[]
for v in vs:
 z=five[v];g=h[v];t=tg[v]
 def add(source,phenotype,beta,se,pval,interpret):
  beta=float(beta);se=float(se) if se not in (None,'','#NA') else None
  out.append({'variant':v,'counted_allele':'G','other_allele':'A','source':source,'phenotype':phenotype,'effect':beta,'se':se if se else '','ci_low':beta-1.95996398454*se if se else '','ci_high':beta+1.95996398454*se if se else '','p_reported':pval,'scale_or_caveat':interpret})
 for label,prefix in [('FI meta-analysis','FI'),('FFS','FFS'),('UKB insomnia','Jansen_insomnia'),('FinnGen insomnia','R13_insomnia')]:
  add(label,label,z[prefix+'_beta'],z[prefix+'_se'],z[prefix+'_p'],'source beta; FFS effect-allele coding inferred' if label=='FFS' else 'source beta')
 add('FinnGen HFRS','HFRS',g['hfrs_beta_on_counted_allele'],g['hfrs_se'],g['hfrs_p'],'HFRS beta, source outcome scale')
 add('UKB FI component','FI',t['UKB_beta_aligned'],t['UKB_se'],t['UKB_p'],'component of FI meta-analysis, not independent')
 add('TwinGene FI component','FI',t['TG_beta_aligned_adapter'],t['TG_se'],t['TG_p'],'component of FI meta-analysis; sign adapter inferred, not author confirmed')
 rr=raw[v]; odds=float(rr['effect_value']);odds=odds if rr['A1']=='G' else 1/odds
 out.append({'variant':v,'counted_allele':'G','other_allele':'A','source':'MVP insomnia','phenotype':'MVP insomnia','effect':odds,'se':'','ci_low':'','ci_high':'','p_reported':rr['p'],'scale_or_caveat':'odds ratio; source rounds OR and does not report SE; CI unavailable'})
q=p/'top_two_candidate_effects.csv'
with q.open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=out[0]);w.writeheader();w.writerows(out)
print('hash',hashlib.sha256(q.read_bytes()).hexdigest())
for r in out:print(r['variant'],r['source'],r['effect'],r['ci_low'],r['ci_high'],r['p_reported'])
