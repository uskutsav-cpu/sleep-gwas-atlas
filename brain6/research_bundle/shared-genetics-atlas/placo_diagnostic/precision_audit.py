import csv,gzip,pathlib,collections,json,math,statistics
root=pathlib.Path('/Users/swethasunilkumar/Documents/Codex/2026-10-27/brain6/outputs/shared-genetics-atlas/placo_diagnostic'); sample={r['SNP']:r for r in csv.DictReader((root/'random_controls.tsv').open(),delimiter='\t')};paths={'S':'/Volumes/Extreme SSD/sleep-frailty-full-autonomous-v5-20261009/frailty_paper/full_autonomous_research_v5/genetics/new_public_validation/filtered/R13_INSOMNIA_HM3.tsv','FI':'/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/data/raw/frailty.txt.gz','FFS':'/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/frailty_paper/data/gwas/fried_frailty_score/Frailty_GWAS_maf0.01_info0.9_hwe1e-6_FUMA.tsv'}
cols={'S':('SNP','A1','A2','BETA','SE','P','Z_insomnia'),'FI':('variant_id','effect_allele','other_allele','beta','standard_error','p_value','Z_FI'),'FFS':('SNP','A1','A2','BETA','SE','P','Z_FFS')}
output=[]; counts={}
for key,path in paths.items():
 seen=set();count=collections.Counter();ks,ka,kb,kbeta,kse,kp,kz=cols[key]
 with (gzip.open(path,'rt') if path.endswith('.gz') else open(path)) as f:
  for r in csv.DictReader(f,delimiter='\t'):
   count['source_rows']+=1;snp=r[ks]
   if snp not in sample:continue
   if snp in seen:count['duplicate_rsid_rows']+=1;continue
   seen.add(snp);m=sample[snp];a=r[ka].upper();b=r[kb].upper();c=m['A1'];d=m['A2']
   if (a,b)==(c,d): sign=1
   elif (a,b)==(d,c):sign=-1
   else:count['allele_mismatch']+=1;continue
   try:
    beta=float(r[kbeta]);se=float(r[kse]);p=float(r[kp]);z=float(m[kz]);assert se>0 and math.isfinite(beta) and math.isfinite(se) and 0<p<=1
    wald=sign*beta/se
    pred=math.copysign(-statistics.NormalDist().inv_cdf(p/2),sign*beta) if p<1 else 0
   except (ValueError,AssertionError,ZeroDivisionError): count['invalid']+=1;continue
   count['matched']+=1
   output.append({'source':key,'SNP':snp,'A1':c,'A2':d,'source_A1':a,'source_A2':b,'source_beta':beta,'source_se':se,'source_p':p,'z_munged':z,'z_wald_beta_se':wald,'z_from_reported_p_and_beta_sign':pred,'delta_munged_minus_wald':z-wald,'delta_munged_minus_p_derived':z-pred})
 counts[key]=dict(count);print(key,dict(count),flush=True)
with (root/'precision_audit_sample.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=output[0],delimiter='\t');w.writeheader();w.writerows(output)
summary={}
for k in paths:
 a=[r for r in output if r['source']==k]; ds=sorted(abs(r['delta_munged_minus_wald']) for r in a);dp=sorted(abs(r['delta_munged_minus_p_derived']) for r in a)
 def quant(x,q):return x[min(len(x)-1,int(q*(len(x)-1)))]
 summary[k]={'counts':counts[k],'median_abs_z_minus_beta_se':quant(ds,.5),'p95_abs_z_minus_beta_se':quant(ds,.95),'max_abs_z_minus_beta_se':max(ds),'n_abs_z_minus_beta_se_gt_0.01':sum(x>.01 for x in ds),'median_abs_z_minus_p_derived':quant(dp,.5),'p95_abs_z_minus_p_derived':quant(dp,.95),'max_abs_z_minus_p_derived':max(dp),'n_abs_z_minus_p_derived_gt_0.001':sum(x>.001 for x in dp)}
(root/'precision_audit_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
