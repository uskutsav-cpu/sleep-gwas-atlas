import csv,gzip,collections,json,hashlib,pathlib,math
from pyliftover import LiftOver
q=list(csv.DictReader(open('work/stage3_QTD000285_ADARB1_splice_1mb.tsv'),delimiter='\t'));qmap=collections.defaultdict(list)
for x in q:
 rs=x['rsid'];pair={x['ref'],x['alt']}
 if rs and rs!='NA' and len(pair)==2 and pair<=set('ACGT') and pair not in ({'A','T'},{'C','G'}) and 0.01<=float(x['maf'])<=0.5 and float(x['se'])>0:qmap[rs].append(x)
qmap={k:v[0] for k,v in qmap.items() if len(v)==1};print('eligible unique qtl rsIDs',len(qmap),flush=True)
root=pathlib.Path('/Volumes/Extreme SSD/Codex Archive/mak2023_same_cohort');base=pathlib.Path('/Volumes/Extreme SSD/sleep-frailty-full-autonomous-v5-20261009/frailty_paper/full_autonomous_research_v5/genetics/new_public_validation/raw');sources={'Mak_FI':root/'fi.fastGWA','Mak_FP':root/'fp.fastGWA','Ye_FFS':pathlib.Path('/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/frailty_paper/data/gwas/fried_frailty_score/Frailty_GWAS_maf0.01_info0.9_hwe1e-6_FUMA.tsv'),'FinnGen_insomnia':base/'finngen_R13_F5_INSOMNIA.gz'};lo=LiftOver('work/hg19ToHg38.over.chain.gz'); results={}; summary={}
for name,path in sources.items():
 found=collections.defaultdict(list); scanned=0
 with (gzip.open(path,'rt') if str(path).endswith('.gz') else path.open()) as f:
  r=csv.DictReader(f,delimiter='\t')
  for row in r:
   scanned+=1;ids=row['rsids'].split(',') if name=='FinnGen_insomnia' else [row['SNP']]
   for rs in ids:
    if rs in qmap:found[rs].append(row)
 print(name,'scanned',scanned,'ids',len(found),flush=True);aligned=[];reasons=collections.Counter()
 for rs,rowlist in found.items():
  if len(rowlist)!=1:reasons['duplicate_gwas_id']+=1;continue
  r=rowlist[0];v=qmap[rs];
  try:
   if name=='FinnGen_insomnia':chrom=r['#chrom'];pos=int(r['pos']);a=r['alt'];b=r['ref'];beta=float(r['beta']);se=float(r['sebeta']);af=float(r['af_alt']);mapped=pos;strand='+'
   else:
    chrom=r['CHROM'] if name=='Ye_FFS' else r['CHR'];pos=int(r['POS']);a=r['A1'];b=r['A2'];beta=float(r['BETA']);se=float(r['SE']);af=None if name=='Ye_FFS' else float(r['AF1']); hits=lo.convert_coordinate('chr'+chrom,pos-1) or [];hits=[x for x in hits if x[0]=='chr'+chrom and x[2]=='+'];
    if len(hits)!=1:reasons['liftover_nonunique_or_negative']+=1;continue
    mapped=hits[0][1]+1;strand='+'
   if chrom!='21' or mapped!=int(v['position']) or {a,b}!={v['ref'],v['alt']}:reasons['coord_or_allele_mismatch']+=1;continue
   if not (math.isfinite(beta) and math.isfinite(se) and se>0 and (af is None or .01<=af<=.99)):reasons['gwas_qc']+=1;continue
   qb=float(v['beta']);qs=float(v['se']);assert math.isfinite(qb) and math.isfinite(qs) and qs>0
   aligned.append({'variant':rs,'position_GRCh38':mapped,'ref':v['ref'],'alt':v['alt'],'qtl_beta_alt':qb,'qtl_se':qs,'qtl_p':float(v['pvalue']),'qtl_maf':float(v['maf']),'gwas_beta_alt':beta if a==v['alt'] else -beta,'gwas_se':se,'gwas_p':float(r['pval'] if name=='FinnGen_insomnia' else r['P']),'gwas_af_alt':(af if a==v['alt'] else 1-af) if af is not None else None})
  except (ValueError,TypeError,KeyError,AssertionError):reasons['invalid_fields']+=1
 aligned.sort(key=lambda x:x['position_GRCh38']);out=pathlib.Path('work/stage3_adarb1_muscle_splice_'+name+'_aligned.csv')
 with out.open('w',newline='') as f:
  fields=['variant','position_GRCh38','ref','alt','qtl_beta_alt','qtl_se','qtl_p','qtl_maf','gwas_beta_alt','gwas_se','gwas_p','gwas_af_alt'];w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(aligned)
 summary[name]={'source':str(path),'scanned':scanned,'aligned':len(aligned),'matched_ids':len(found),'exclusions':dict(reasons),'aligned_sha256':hashlib.sha256(out.read_bytes()).hexdigest(),'candidate_present':any(x['variant']=='rs2838788' for x in aligned)}
 print(name,summary[name],flush=True)
pathlib.Path('outputs/stage3_adarb1_muscle_splice_harmonization_receipt.json').write_text(json.dumps(summary,indent=2)+'\n')
