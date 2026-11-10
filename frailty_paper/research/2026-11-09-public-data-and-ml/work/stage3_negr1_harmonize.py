import csv,gzip,collections,json,hashlib,pathlib,math
from pyliftover import LiftOver
w=[x for x in csv.DictReader(open('work/stage3_negr1_worklist.csv')) if x['candidate']];load=json.load(open('work/stage3_negr1_extract_summary.json'));ok={x['candidate'] for x in load['items'] if x['status']=='extracted' and x['gene_rows']>0};w=[x for x in w if x['candidate'] in ok]
qmaps={}; reverse=collections.defaultdict(list)
for x in w:
 c=x['candidate'];qmap=collections.defaultdict(list);f=pathlib.Path('work/stage3_negr1_qtl')/(c+'.tsv')
 for v in csv.DictReader(f.open(),delimiter='\t'):
  try:
   pair={v['ref'],v['alt']};good=v['rsid'] and v['rsid']!='NA' and len(pair)==2 and pair<=set('ACGT') and pair not in ({'A','T'},{'C','G'}) and .01<=float(v['maf'])<=.5 and float(v['se'])>0
  except ValueError:good=False
  if good:qmap[v['rsid']].append(v)
 qmap={k:v[0] for k,v in qmap.items() if len(v)==1};qmaps[c]=qmap
 for rs in qmap:reverse[rs].append(c)
print('tasks',len(w),'unique ids',len(reverse),flush=True);lo=LiftOver('work/hg19ToHg38.over.chain.gz');traits=('Atkins_FI','Jansen_insomnia','Mak_FI','Mak_FP','Ye_FFS','FinnGen_insomnia');summary={};outdir=pathlib.Path('work/stage3_negr1_qtl/aligned');outdir.mkdir(exist_ok=True)
for name in traits:
 path=pathlib.Path('work/stage3_region_cache')/(name+'.regional.tsv.gz');found=collections.defaultdict(list);scanned=0
 with gzip.open(path,'rt') as f:
  for row in csv.DictReader(f,delimiter='\t'):
   scanned+=1;rsids=row['rsids'].split(',') if name=='FinnGen_insomnia' else [row['SNP']]
   for rs in rsids:
    if rs in reverse:
     for c in reverse[rs]:found[(c,rs)].append(row)
 for item in w:
  c=item['candidate'];chrom=item['marker_grch38'].split(':')[0].removeprefix('chr');aligned=[];reasons=collections.Counter();matches=0
  for rs,v in qmaps[c].items():
   rr=found.get((c,rs),[])
   if not rr:continue
   matches+=1
   if len(rr)!=1:reasons['duplicate_gwas_id']+=1;continue
   r=rr[0]
   try:
    if name=='FinnGen_insomnia':ch=r['#chrom'];pos=int(r['pos']);a=r['alt'];b=r['ref'];beta=float(r['beta']);se=float(r['sebeta']);af=float(r['af_alt']);mapped=pos;p=float(r['pval']);N=''
    else:
     ch=r['CHROM'] if name=='Ye_FFS' else r['CHR'];pos=int(r['BP'] if name in ('Atkins_FI','Jansen_insomnia') else r['POS']);a=r['A1'];b=r['A2'];beta=float(r['BETA']);se=float(r['SE']);af=None if name=='Ye_FFS' else float(r['FRQ'] if name in ('Atkins_FI','Jansen_insomnia') else r['AF1']);p=float(r['P']);N=r.get('N','');hits=lo.convert_coordinate('chr'+ch,pos-1) or [];hits=[z for z in hits if z[0]=='chr'+ch and z[2]=='+']
     if len(hits)!=1:reasons['liftover_nonunique_or_negative']+=1;continue
     mapped=hits[0][1]+1
    if ch!=chrom or mapped!=int(v['position']) or {a,b}!={v['ref'],v['alt']}:reasons['coordinate_or_allele_mismatch']+=1;continue
    if not(math.isfinite(beta) and math.isfinite(se) and se>0 and (af is None or .01<=af<=.99)):reasons['gwas_qc']+=1;continue
    qb=float(v['beta']);qs=float(v['se']);assert math.isfinite(qb) and math.isfinite(qs) and qs>0
    aligned.append({'variant':rs,'position_GRCh38':mapped,'ref':v['ref'],'alt':v['alt'],'qtl_beta_alt':qb,'qtl_se':qs,'qtl_p':float(v['pvalue']),'qtl_maf':float(v['maf']),'gwas_beta_alt':beta if a==v['alt'] else -beta,'gwas_se':se,'gwas_p':p,'gwas_af_alt':(af if a==v['alt'] else 1-af) if af is not None else None,'gwas_N':N})
   except (ValueError,TypeError,KeyError,AssertionError):reasons['invalid_fields']+=1
  aligned.sort(key=lambda z:z['position_GRCh38']);out=outdir/(c+'_'+name+'.csv');fields=['variant','position_GRCh38','ref','alt','qtl_beta_alt','qtl_se','qtl_p','qtl_maf','gwas_beta_alt','gwas_se','gwas_p','gwas_af_alt','gwas_N']
  with out.open('w',newline='') as f:writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader();writer.writerows(aligned)
  summary[c+'|'+name]={'n':len(aligned),'matched_ids':matches,'candidate_present':any(z['variant']==c for z in aligned),'exclusions':dict(reasons),'output_sha256':hashlib.sha256(out.read_bytes()).hexdigest()}
 print(name,'cache rows',scanned,flush=True)
pathlib.Path('outputs/stage3_negr1_qtl_harmonization_receipt.json').write_text(json.dumps({'results':summary,'caveat':'regional caches; unique rsID, nonpalindromic biallelic alleles, QTL MAF >=1%, coordinate liftover and GWAS allele pair matching; input builds, ancestry and quality not independently established'},indent=2)+'\n')
