import csv,math,json,hashlib,pathlib
p=pathlib.Path('outputs')
def load(f,k):return {x[k]:x for x in csv.DictReader((p/f).open())}
j=load('song124_mvp_hfrs_joint_screen.csv','variant');meta=load('published_candidate_marginal_lookup.csv','candidate');genes=load('song124_22_nominal_positional_genes.csv','variant'); old=load('unvalidated_15_variant_score_prototype.csv','variant');ld=json.load(open(p/'song124_22_nominal_reference_ld.json'));out=[]
for v,x in j.items():
 if not x.get('joint_p_max') or float(x['joint_p_max'])>=.05:continue
 y=meta[v];g=genes[v];a=x['HFRS_effect_allele'];b=x['HFRS_other_allele'];d={'variant':v,'GRCh37_chr_assumed':x['song_chr'],'GRCh37_pos_assumed':x['song_position_reported'],'GRCh38_chr':x['GRCh38_chr'],'GRCh38_pos':x['GRCh38_pos'],'counted_allele':a,'other_allele':b,'nearest_gene_body_positional':g.get('nearest_names',''),'nearest_gene_biotype':g.get('nearest_biotypes',''),'distance_bp':g.get('distance_bp',''),'reference_alleles':g.get('reference_alleles',''),'1000G_EUR_matching_biallelic_contrast':ld['variant_status'][v]['status'],'in_original_15':v in old}
 for source,prefix in [('FI','FI'),('FFS','FFS'),('UKB_insomnia','insomnia')]:
  aa,bb=y[prefix+'_A1'],y[prefix+'_A2'];assert {aa,bb}=={a,b},(v,source,aa,bb,a,b); assert y[prefix+'_same_reported_position_as_song']=='True'
  d[source+'_beta_aligned']=float(y[prefix+'_BETA'])*(1 if aa==a else -1);d[source+'_p']=float(y[prefix+'_P'])
 d['HFRS_beta_aligned']=float(x['HFRS_effect']);d['HFRS_se']=float(x['HFRS_se']);d['HFRS_p']=float(x['HFRS_p']);d['MVP_insomnia_log_OR_aligned_from_rounded_OR']=math.log(float(x['MVP_effect']))*(1 if x['MVP_effect_allele']==a else -1);d['MVP_insomnia_p']=float(x['MVP_p']);d['joint_max_p']=float(x['joint_p_max']);d['numerical_bonferroni_124']=float(x['joint_bonferroni_124']);d['all_five_nominal']=all(d[k]<.05 for k in ('FI_p','FFS_p','UKB_insomnia_p','HFRS_p','MVP_insomnia_p'));d['all_five_same_sign']=len({math.copysign(1,d[k]) for k in ('FI_beta_aligned','FFS_beta_aligned','UKB_insomnia_beta_aligned','HFRS_beta_aligned','MVP_insomnia_log_OR_aligned_from_rounded_OR')})==1
 out.append(d)
out.sort(key=lambda x:x['joint_max_p']);assert len(out)==22
q=p/'expanded_22_candidate_risk_gene_table.csv'
with q.open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=out[0]);w.writeheader();w.writerows(out)
print('sha256',hashlib.sha256(q.read_bytes()).hexdigest(),'five nominal',sum(x['all_five_nominal'] for x in out),'five signs',sum(x['all_five_same_sign'] for x in out),'both',sum(x['all_five_nominal'] and x['all_five_same_sign'] for x in out))
print('multiallelic',sum(len(x['reference_alleles'].split('/'))>2 for x in out))
for x in out:print(x['variant'],x['nearest_gene_body_positional'],x['counted_allele'],x['FI_p'],x['FFS_p'],x['all_five_same_sign'])
