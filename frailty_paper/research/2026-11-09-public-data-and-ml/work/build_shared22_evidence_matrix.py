import csv,pathlib,json,math,hashlib
load=lambda p:{r['variant']:r for r in csv.DictReader(open(p))}
a=load('outputs/expanded_22_candidate_risk_gene_table.csv');b=load('outputs/mak2023_22_candidate_lookup.csv');c=load('outputs/22_finngen_insomnia_followup.csv');score=list(csv.DictReader(open('outputs/unvalidated_22_shared_direction_score.csv')));assert len(score)==len(a)==len(b)==len(c)==22
spec=[('atkins_FI','FI_beta_aligned','FI_p','a'),('ye_FFS','FFS_beta_aligned','FFS_p','a'),('mak_FI','mak_FI_beta','mak_FI_p','b'),('mak_FP','mak_FP_beta','mak_FP_p','b'),('hfrs','HFRS_beta_aligned','HFRS_p','a'),('ukb_insomnia','UKB_insomnia_beta_aligned','UKB_insomnia_p','a'),('mvp_insomnia','MVP_insomnia_log_OR_aligned_from_rounded_OR','MVP_insomnia_p','a'),('finngen_insomnia','r13_beta','r13_p','c')];out=[]
for r in score:
 rs=r['variant'];a0,b0,c0=a[rs],b[rs],c[rs];assert all(x['counted_allele']==a0['counted_allele'] and x['other_allele']==a0['other_allele'] for x in (b0,c0));assert {r['counted_allele'],r['other_allele']}=={a0['counted_allele'],a0['other_allele']};sign=1 if r['counted_allele']==a0['counted_allele'] else -1
 d={'variant':rs,'chr_GRCh37_assumed':r['chr_GRCh37_assumed'],'position_GRCh37_assumed':r['position_GRCh37_assumed'],'chr_GRCh38':a0['GRCh38_chr'],'position_GRCh38':a0['GRCh38_pos'],'score_counted_allele':r['counted_allele'],'other_allele':r['other_allele'],'score_weight':1,'nearest_gene_body_positional':r['nearest_gene_body_positional'],'reference_alleles':r['reference_alleles']}
 n=0;positive=0
 for name,beta,p,table in spec:
  x={'a':a0,'b':b0,'c':c0}[table];v=float(x[beta])*sign;pv=float(x[p]);assert math.isfinite(v) and 0<=pv<=1
  d[name+'_beta_aligned']=v;d[name+'_p']=pv;n+=pv<.05;positive+=v>0
 for name,field,table in [('mak_FI','mak_FI_se','b'),('mak_FP','mak_FP_se','b'),('hfrs','HFRS_se','a'),('finngen_insomnia','r13_se','c')]:
  x={'a':a0,'b':b0,'c':c0}[table];se=float(x[field]);v=d[name+'_beta_aligned'];d[name+'_se']=se;d[name+'_ci95_low']=v-1.95996398454*se;d[name+'_ci95_high']=v+1.95996398454*se
 d.update(n_sources_nominal=n,n_sources_positive_point=positive,all_eight_nominal_same_sign=(n==8 and positive==8),caveat='Postselected; overlapping or potentially overlapping cohorts; unadjusted p; nearby gene is positional; score unvalidated')
 out.append(d)
p=pathlib.Path('outputs/shared22_eight_source_evidence_matrix.csv')
with p.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=out[0]);w.writeheader();w.writerows(out)
ret={'candidate_count':len(out),'sources':[x[0] for x in spec],'all_eight_nominal_same_sign':[x['variant'] for x in out if x['all_eight_nominal_same_sign']],'n_each_p_nominal':{name:sum(x[name+'_p']<.05 for x in out) for name,_,_,_ in spec},'source_units':{'atkins_FI':'quantile-normalized FI','ye_FFS':'FFS count','mak_FI':'FI percent units','mak_FP':'FP count','hfrs':'hospital frailty score','ukb_insomnia':'source study phenotype units','mvp_insomnia':'log OR reconstructed from rounded source OR','finngen_insomnia':'log OR'},'output_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'warning':'Descriptive post-selection and source-overlap subject to allele/build and model caveats, not eight independent replications or risk validation.'};pathlib.Path('outputs/shared22_eight_source_evidence_matrix.json').write_text(json.dumps(ret,indent=2)+'\n');print(json.dumps(ret,indent=2))
