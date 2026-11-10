import csv,math,pathlib,json,statistics,ast,gzip,hashlib
p=pathlib.Path;ns={'math':math};tree=ast.parse(p('work/stage3_smim32_coloc_abf.py').read_text());exec(compile(ast.Module(body=[x for x in tree.body if isinstance(x,ast.FunctionDef)],type_ignores=[]),'functions','exec'),ns)
q=[x for x in csv.DictReader(open('work/stage3_QTD000301_FOXP2_2mb.tsv'),delimiter='\t') if x['variant']=='chr7_114534306_T_TA'];assert len(q)==2 and len({(x['beta'],x['se'],x['pvalue'],x['maf']) for x in q})==1;q=q[0]
indels={(x['variant_GRCh38'],x['trait_source']):x for x in csv.DictReader(open('outputs/stage5_indel_grch37_gwas_matches.csv'))};f=next(x for x in csv.DictReader(open('outputs/stage5_finngen_insomnia_exact_queue.csv')) if x['variant_GRCh38']=='chr7:114534306:T>TA');old=list(csv.DictReader(open('outputs/stage3_foxp2_pituitary_coloc_sensitivity.csv')));cs={x['variant'] for x in csv.DictReader(open('outputs/stage3_foxp2_pituitary_susie_credible_set.csv'))};rows=[];coverage=[]
for trait in ('Mak_FI','Mak_FP','Ye_FFS','FinnGen_insomnia'):
 src=p('work/stage3_foxp2_pituitary_'+trait+'_aligned.csv');d=list(csv.DictReader(src.open()));assert not any(x['position_GRCh38']=='114534306' and x['ref']=='T' and x['alt']=='TA' for x in d)
 oldrow=next(x for x in old if x['trait']==trait and float(x['p12'])==1e-5 and float(x['prior_scale'])==1);sd=float(oldrow['gwas_prior_sd'])
 before=ns['posterior']([ns['bf'](float(x['qtl_beta_alt']),float(x['qtl_se']),.15) for x in d],[ns['bf'](float(x['gwas_beta_alt']),float(x['gwas_se']),sd) for x in d],1e-4,1e-4,1e-5);assert max(abs(before[i]-float(oldrow['H'+str(i)])) for i in range(5))<1e-10
 if trait=='FinnGen_insomnia':assert f['status']=='exact REF/ALT match';beta=f['beta_alt'];se=f['se'];pv=f['p_value'];af=f['af_alt']
 else:
  z=indels[('chr7:114534306:T>TA',trait)];assert z['match_status']=='unique_exact_allele_match';beta=z['beta_alt'];se=z['se'];pv=z['p'];af=''
  if trait!='Ye_FFS':
   with gzip.open('work/stage3_region_cache/'+trait+'.regional.tsv.gz','rt') as o:
    hits=[x for x in csv.DictReader(o,delimiter='\t') if x['CHR']=='7' and x['POS']=='114174361' and {x['A1'],x['A2']}=={'T','TA'}]
   assert len(hits)==1;x=hits[0];af=float(x['AF1']) if x['A1']=='TA' else 1-float(x['AF1']);assert .01<float(af)<.99
 v={'variant':'rs66823671','position_GRCh38':'114534306','ref':'T','alt':'TA','qtl_beta_alt':q['beta'],'qtl_se':q['se'],'qtl_p':q['pvalue'],'qtl_maf':q['maf'],'gwas_beta_alt':beta,'gwas_se':se,'gwas_p':pv,'gwas_af_alt':af};assert .01<float(q['maf'])<.5
 combined=d+[v];fn=p('outputs/stage5_foxp2_pituitary_'+trait+'_with_indel_aligned.csv')
 with fn.open('w',newline='') as o:w=csv.DictWriter(o,fieldnames=list(d[0]));w.writeheader();w.writerows(combined)
 support={'chr7_'+x['position_GRCh38']+'_'+x['ref']+'_'+x['alt'] for x in combined};coverage.append({'trait':trait,'previous_n':len(d),'updated_n':len(combined),'credible_set_total':len(cs),'credible_set_matched_before':len(cs&support)-1,'credible_set_matched_after':len(cs&support),'indel_qtl_beta_alt':q['beta'],'indel_qtl_se':q['se'],'indel_gwas_beta_alt':beta,'indel_gwas_se':se,'indel_gwas_p':pv,'gwas_af_alt':af,'source_aligned_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'updated_aligned_sha256':hashlib.sha256(fn.read_bytes()).hexdigest()})
 for scale in (.5,1.,2.):
  l1=[ns['bf'](float(x['qtl_beta_alt']),float(x['qtl_se']),.15) for x in combined];l2=[ns['bf'](float(x['gwas_beta_alt']),float(x['gwas_se']),sd*scale) for x in combined]
  for p12 in (1e-6,1e-5,5e-5):
   h=ns['posterior'](l1,l2,1e-4,1e-4,p12);rows.append({'trait':trait,'n':len(combined),'gwas_prior_sd':sd*scale,'prior_scale':scale,'p12':p12,**{'H'+str(i):h[i] for i in range(5)}})
with open('outputs/stage5_foxp2_indel_coloc_sensitivity.csv','w',newline='') as o:w=csv.DictWriter(o,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
with open('outputs/stage5_foxp2_indel_coloc_coverage.csv','w',newline='') as o:w=csv.DictWriter(o,fieldnames=list(coverage[0]));w.writeheader();w.writerows(coverage)
print('baseline with insertion:');[print(x['trait'],x['n'],'H3',round(x['H3'],6),'H4',round(x['H4'],6)) for x in rows if x['prior_scale']==1 and x['p12']==1e-5];print('coverage',[(x['trait'],x['credible_set_matched_before'],x['credible_set_matched_after']) for x in coverage])
