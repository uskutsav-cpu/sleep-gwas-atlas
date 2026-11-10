import csv,math,collections,json,pathlib,hashlib,runpy
p=pathlib.Path('work/stage3_QTD000285_SMAD5_splice_1mb.tsv');raw=list(csv.DictReader(p.open(),delimiter='\t'));d=collections.defaultdict(list)
for x in raw:
 pair={x['ref'],x['alt']}
 if x['rsid']!='NA' and len(pair)==2 and pair<=set('ACGT') and pair not in ({'A','T'},{'C','G'}) and 0.01<=float(x['maf'])<=.5 and float(x['se'])>0:d[x['rsid']].append(x)
q={k:v[0] for k,v in d.items() if len(v)==1}; f=runpy.run_path('work/stage3_smim32_coloc_abf.py'); bf=f['bf'];post=f['posterior'];diags=json.load(open('outputs/stage3_smim32_mak_coloc_receipt.json'))['inputs'];out=[];rec={}
for name in ('Mak_FI','Mak_FP','FinnGen_insomnia'):
 source=pathlib.Path('work/stage3_smim32_'+name+('_with_n.csv' if name.startswith('Mak') else '_aligned.csv'));r=list(csv.DictReader(source.open()));merged=[]
 for x in r:
  y=q.get(x['variant']);
  if y is None or x['position_GRCh38']!=y['position'] or x['ref']!=y['ref'] or x['alt']!=y['alt']:continue
  merged.append({**x,'qtl_beta_alt':y['beta'],'qtl_se':y['se'],'qtl_p':y['pvalue'],'qtl_maf':y['maf']})
 assert len(merged)>100 and len({x['variant'] for x in merged})==len(merged)
 p_out=pathlib.Path('work/stage3_smad5_'+name+'_aligned.csv')
 with p_out.open('w',newline='') as f2:w=csv.DictWriter(f2,fieldnames=merged[0].keys());w.writeheader();w.writerows(merged)
 rec[name]={'variants':len(merged),'candidate_present':any(x['variant']=='rs7734045' for x in merged),'gwas_prior_sd':.2 if name=='FinnGen_insomnia' else .15*diags[name]['estimated_gwas_sdY'],'matched_sha256':hashlib.sha256(p_out.read_bytes()).hexdigest()}
 for p12 in (1e-6,1e-5,5e-5):
  for factor in (.5,1,2):
   qprior=.15*factor;gp=rec[name]['gwas_prior_sd'];l1=[bf(float(x['qtl_beta_alt']),float(x['qtl_se']),qprior) for x in merged];l2=[bf(float(x['gwas_beta_alt']),float(x['gwas_se']),gp) for x in merged];h=post(l1,l2,1e-4,1e-4,p12)
   out.append({'trait':name,'n_variants':len(merged),'qtl_prior_effect_sd':qprior,'gwas_prior_sd':gp,'p1':1e-4,'p2':1e-4,'p12':p12,'H0':h[0],'H1':h[1],'H2':h[2],'H3':h[3],'H4':h[4],'H3_plus_H4':h[3]+h[4],'H4_cond_H3_H4':h[4]/(h[3]+h[4])})
outfile=pathlib.Path('outputs/stage3_smad5_splice_coloc_abf_sensitivity.csv')
with outfile.open('w',newline='') as f2:w=csv.DictWriter(f2,fieldnames=out[0].keys());w.writeheader();w.writerows(out)
pathlib.Path('outputs/stage3_smad5_splice_coloc_receipt.json').write_text(json.dumps({'event':'5:136132962:136153592:clu_31306_+','qtl_source':'eQTL Catalogue r7 QTD000285 cc filtered; selected event','traits':rec,'output_sha256':hashlib.sha256(outfile.read_bytes()).hexdigest(),'limits':'postselected splice event, inherited eQTL-GWAS intersection restricts support, single-signal custom ABF, missing SNPs, assumptions and sample overlap not fully verified'},indent=2)+'\n')
for x in out:
 if x['qtl_prior_effect_sd']==.15 and x['p12']==1e-5:print(x)
