import csv,math,pathlib,hashlib
p=pathlib.Path('outputs')
def get(file,key):return {r[key]:r for r in csv.DictReader((p/file).open())}
a=get('song124_mvp_hfrs_joint_screen.csv','variant');b=get('published_candidate_marginal_lookup.csv','candidate');c=get('six_joint_hits_finngen_insomnia_followup.csv','variant');g=get('song124_six_hit_positional_genes.csv','variant')
out=[]
for variant in c:
 x=a[variant];y=b[variant];z=c[variant];counted=x['HFRS_effect_allele'];other=x['HFRS_other_allele'];o={'variant':variant,'counted_allele':counted,'other_allele':other,'positional_gene':g[variant].get('nearest_names',''),'reference_alleles':g[variant].get('reference_alleles','')}
 for prefix in ('FI','FFS','insomnia'):
  aa=y.get(prefix+'_A1','');bb=y.get(prefix+'_A2','');bet=y.get(prefix+'_BETA','');ps=y.get(prefix+'_P','')
  if aa and bb and bet and {aa,bb}=={counted,other}:
   o[prefix+'_aligned_beta']=float(bet)*(1 if aa==counted else -1);o[prefix+'_p']=ps
  else:o[prefix+'_aligned_beta']='';o[prefix+'_p']=ps
 o['HFRS_beta']=float(x['HFRS_effect']);o['HFRS_se']=float(x['HFRS_se']);o['HFRS_p']=x['HFRS_p']
 odds=float(x['MVP_effect']);o['MVP_logOR_aligned']=math.log(odds)*(1 if x['MVP_effect_allele']==counted else -1);o['MVP_p']=x['MVP_p']
 o['R13_beta_aligned']=z.get('r13_beta_aligned','');o['R13_se']=z.get('r13_se','');o['R13_p']=z.get('r13_p','')
 vals=[o.get(k,'') for k in ('FI_aligned_beta','FFS_aligned_beta','insomnia_aligned_beta','HFRS_beta','MVP_logOR_aligned','R13_beta_aligned')]
 o['all_six_positive_or_all_six_negative']=all(v!='' for v in vals) and (all(float(v)>0 for v in vals) or all(float(v)<0 for v in vals));o['all_six_nominal_p_lt_005']=all(o.get(k,'') and float(o[k])<.05 for k in ('FI_p','FFS_p','insomnia_p','HFRS_p','MVP_p','R13_p'));out.append(o)
q=p/'six_candidates_six_trait_direction_matrix.csv';keys=list(dict.fromkeys(k for r in out for k in r))
with q.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(out)
print('sha256',hashlib.sha256(q.read_bytes()).hexdigest())
for x in out:print(x['variant'],x['positional_gene'],'direction_same',x['all_six_positive_or_all_six_negative'],'all_p',x['all_six_nominal_p_lt_005'],'r13_p',x['R13_p'],'ffs_p',x['FFS_p'])
