import csv,gzip,math,json,hashlib,ast
from pathlib import Path
from pyliftover import LiftOver
root=Path(__file__).resolve().parents[1];ns={'math':math};t=ast.parse((root/'work/stage3_smim32_coloc_abf.py').read_text());exec(compile(ast.Module(body=[n for n in t.body if isinstance(n,ast.FunctionDef)],type_ignores=[]),'coloc','exec'),ns)
q=list(csv.DictReader((root/'work/tenk10k_cd4_smad5_peak_full_cis.csv').open()));lo=LiftOver(str(root/'work/hg38ToHg19.over.chain.gz'));eligible={};stats={'all_qtl':len(q),'non_snv':0,'low_af':0,'liftover_nonunique':0,'eligible':0}
for r in q:
 c,pos,ref,alt=r['variant_id'].split(':');pos=int(pos); af=float(r['af'])
 if len(ref)!=1 or len(alt)!=1:stats['non_snv']+=1;continue
 if not 0.01<=af<=0.99:stats['low_af']+=1;continue
 m=lo.convert_coordinate('chr'+c,pos-1)
 if len(m)!=1 or m[0][0]!='chr5' or m[0][2]!='+':stats['liftover_nonunique']+=1;continue
 key=(int(m[0][1])+1,frozenset((ref,alt)));eligible[key]=r
stats['eligible']=len(eligible); print('eligible',stats,flush=True)
results=[];coverage=[]
for trait,sd in [('Mak_FI',1.1004195728364333),('Mak_FP',.1232461359076205),('FinnGen_insomnia',.2)]:
 path=root/f'work/stage3_region_cache/{trait}.regional.tsv.gz'; matches={};dups=set()
 with gzip.open(path,'rt') as f:
  for x in csv.DictReader(f,delimiter='\t'):
   if trait=='FinnGen_insomnia':
    if x['#chrom']!='5':continue
    bpos=int(x['pos']);ref=x['ref'];alt=x['alt']; map37=lo.convert_coordinate('chr5',bpos-1)
    if len(map37)!=1 or map37[0][0]!='chr5' or map37[0][2]!='+':continue
    key=(int(map37[0][1])+1,frozenset((ref,alt)));effect=alt;af=x['af_alt'];beta=x['beta'];se=x['sebeta'];gpos=bpos
   else:
    if x['CHR']!='5':continue
    key=(int(x['POS']),frozenset((x['A1'],x['A2'])));effect=x['A1'];af=x['AF1'];beta=x['BETA'];se=x['SE'];gpos=None
   if key not in eligible:continue
   z=eligible[key];c,pos,qref,qalt=z['variant_id'].split(':')
   if trait=='FinnGen_insomnia' and (gpos!=int(pos) or ref!=qref or alt!=qalt):continue
   if not all([af,beta,se]) or any(v in ('NA','nan','.') for v in [af,beta,se]):continue
   freq=float(af) if effect==qalt else 1-float(af)
   if not .01<=freq<=.99 or float(se)<=0 or abs(freq-float(z['af']))>.2:continue
   if {qref,qalt} in ({'A','T'},{'C','G'}) and (0.4 <= float(z['af']) <= 0.6 or abs(freq-float(z['af']))>.1):continue
   entry={'variant_GRCh38':z['variant_id'],'gwas_pos37':key[0],'qtl_alt_af':z['af'],'qtl_beta_alt':z['slope'],'qtl_se':z['slope_se'],'qtl_p':z['pval_nominal'],'gwas_af_alt':freq,'gwas_beta_alt':float(beta)*(1 if effect==qalt else -1),'gwas_se':se}
   if key in matches:dups.add(key)
   else:matches[key]=entry
 for k in dups:matches.pop(k,None)
 rows=list(matches.values()); assert rows
 fp=root/f'outputs/stage6_tenk10k_peak_pal_freq_{trait}_aligned.csv'
 with fp.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
 cov={'trait':trait,'matched_n':len(rows),'duplicates_excluded':len(dups),'qtl_lead_present':any(x['variant_GRCh38']=='5:136187682:T:G' for x in rows),'candidate_present':any(x['variant_GRCh38']=='5:136202653:A:C' for x in rows),'min_af_delta_median':sorted(abs(float(x['qtl_alt_af'])-float(x['gwas_af_alt'])) for x in rows)[len(rows)//2],'max_abs_frequency_difference':max(abs(float(x['qtl_alt_af'])-float(x['gwas_af_alt'])) for x in rows),'file_sha256':hashlib.sha256(fp.read_bytes()).hexdigest()};coverage.append(cov)
 for qsd in [.075,.15,.30]:
  for p12 in [1e-6,1e-5,5e-5]:
   h=ns['posterior']([ns['bf'](float(x['qtl_beta_alt']),float(x['qtl_se']),qsd) for x in rows],[ns['bf'](float(x['gwas_beta_alt']),float(x['gwas_se']),sd) for x in rows],1e-4,1e-4,p12)
   results.append({'trait':trait,'n':len(rows),'qtl_prior_sd':qsd,'gwas_prior_sd':sd,'p12':p12,**{'H'+str(i):h[i] for i in range(5)}})
 print(trait,cov,flush=True)
out=root/'outputs/stage6_tenk10k_peak_pal_freq_disease_coloc_sensitivity.csv'
with out.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=results[0].keys());w.writeheader();w.writerows(results)
(root/'outputs/stage6_tenk10k_peak_pal_freq_disease_coloc_receipt.json').write_text(json.dumps({'qtl_filter':stats,'coverage':coverage,'limitations':'Exploratory postselected single-signal ABF; regional cache coverage, no QTL SNP-level imputation info, ancestry mismatch, overlap, sample sizes, phenotype scales and effect priors require review. Not causal or independent score validation.'},indent=2)+'\n')
print('baseline',[(r['trait'],r['H3'],r['H4']) for r in results if r['qtl_prior_sd']==.15 and r['p12']==1e-5])
