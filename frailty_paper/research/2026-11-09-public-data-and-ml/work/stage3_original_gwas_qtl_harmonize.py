import csv,gzip,collections,json,hashlib,pathlib,math,sys
from pyliftover import LiftOver
analysis,qtlpath,chrom,candidate=sys.argv[1:];q=list(csv.DictReader(open(qtlpath),delimiter='\t'));qmap=collections.defaultdict(list)
for x in q:
 rs=x['rsid'];pair={x['ref'],x['alt']}
 try:valid=rs and rs!='NA' and len(pair)==2 and pair<=set('ACGT') and pair not in ({'A','T'},{'C','G'}) and .01<=float(x['maf'])<=.5 and float(x['se'])>0
 except ValueError:valid=False
 if valid:qmap[rs].append(x)
qmap={k:v[0] for k,v in qmap.items() if len(v)==1};lo=LiftOver('work/hg19ToHg38.over.chain.gz');summary={}
for name in ('Atkins_FI','Jansen_insomnia'):
 path=pathlib.Path('work/stage3_region_cache')/(name+'.regional.tsv.gz');found=collections.defaultdict(list);scanned=0
 with gzip.open(path,'rt') as f:
  for row in csv.DictReader(f,delimiter='\t'):
   scanned+=1
   if row['SNP'] in qmap:found[row['SNP']].append(row)
 aligned=[];reasons=collections.Counter()
 for rs,rsrows in found.items():
  if len(rsrows)!=1:reasons['duplicate_gwas_id']+=1;continue
  r=rsrows[0];v=qmap[rs]
  try:
   ch=r['CHR'];pos=int(r['BP']);a=r['A1'];b=r['A2'];beta=float(r['BETA']);se=float(r['SE']);af=float(r['FRQ']);hits=lo.convert_coordinate('chr'+ch,pos-1) or [];hits=[z for z in hits if z[0]=='chr'+ch and z[2]=='+']
   if len(hits)!=1:reasons['liftover_nonunique_or_negative']+=1;continue
   mapped=hits[0][1]+1
   if ch!=chrom or mapped!=int(v['position']) or {a,b}!={v['ref'],v['alt']}:reasons['coordinate_or_allele_mismatch']+=1;continue
   if not (math.isfinite(beta) and math.isfinite(se) and se>0 and .01<=af<=.99):reasons['gwas_qc']+=1;continue
   qb=float(v['beta']);qs=float(v['se']);assert math.isfinite(qb) and math.isfinite(qs) and qs>0
   aligned.append({'variant':rs,'position_GRCh38':mapped,'ref':v['ref'],'alt':v['alt'],'qtl_beta_alt':qb,'qtl_se':qs,'qtl_p':float(v['pvalue']),'qtl_maf':float(v['maf']),'gwas_beta_alt':beta if a==v['alt'] else -beta,'gwas_se':se,'gwas_p':float(r['P']),'gwas_af_alt':af if a==v['alt'] else 1-af,'gwas_N':r['N']})
  except (ValueError,TypeError,KeyError,AssertionError):reasons['invalid_fields']+=1
 aligned.sort(key=lambda x:x['position_GRCh38']);out=pathlib.Path('work/stage3_'+analysis+'_'+name+'_aligned.csv')
 with out.open('w',newline='') as f:
  fields=['variant','position_GRCh38','ref','alt','qtl_beta_alt','qtl_se','qtl_p','qtl_maf','gwas_beta_alt','gwas_se','gwas_p','gwas_af_alt','gwas_N'];w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(aligned)
 summary[name]={'cache':str(path),'cache_rows_scanned':scanned,'matching_ids':len(found),'aligned':len(aligned),'candidate_present':any(x['variant'] in candidate.split(',') for x in aligned),'exclusions':dict(reasons),'aligned_sha256':hashlib.sha256(out.read_bytes()).hexdigest()}
 print(analysis,name,summary[name],flush=True)
pathlib.Path('outputs/stage3_'+analysis+'_original_gwas_harmonization_receipt.json').write_text(json.dumps({'qtl_input':qtlpath,'chromosome':chrom,'qtl_eligible_unique_rsids':len(qmap),'results':summary,'caveat':'candidate region cache is source-filtered; original project harmonized files exclude MHC, palindromic variants, MAF<=1%; genome builds and allele semantics not independently authenticated'},indent=2)+'\n')
