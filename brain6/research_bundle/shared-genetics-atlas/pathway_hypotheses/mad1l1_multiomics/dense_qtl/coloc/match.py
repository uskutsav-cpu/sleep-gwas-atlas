import pathlib,csv,collections,math,json,hashlib
p=pathlib.Path(__file__).resolve().parent;qtl=p.parent
spec={'AC110781_3':qtl/'Brain_Nucleus_accumbens_basal_ganglia_ENSG00000176349.tsv','MRM2':qtl/'Brain_Caudate_basal_ganglia_ENSG00000122687.tsv'}
def rows(path):
 with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def snv(a,b):return len(a)==len(b)==1 and a in 'ACGT' and b in 'ACGT' and a!=b and frozenset((a,b)) not in (frozenset('AT'),frozenset('CG'))
def valid(b,s,p,f):return all(math.isfinite(x) for x in (b,s,p,f)) and s>0 and 0<=p<=1 and 0<f<1
outsummary={}
for tissue,path in spec.items():
 q=rows(path);qindex=collections.defaultdict(list);qc=collections.Counter()
 for x in q:
  try:
   ref,alt=x['ref'],x['alt'];beta,se,pv=map(float,(x['beta'],x['se'],x['pvalue']));an,ac=map(float,(x['an'],x['ac']));f=float(x['maf'])
   if not snv(ref,alt) or not valid(beta,se,pv,f) or an<=0:raise ValueError()
  except (ValueError,ZeroDivisionError):qc['qtl_ineligible']+=1;continue
  qc['qtl_eligible_rows']+=1;qindex[(x['chromosome'],x['position'],frozenset((ref,alt)))].append(x)
 for trait in ('insomnia','hfrs'):
  c=qc.copy();gw=rows(p/(trait+'_region.tsv'));index=collections.defaultdict(list)
  for x in gw:
   try:
    if trait=='insomnia':ch,pos,a,b=x['#chrom'],x['pos'],x['alt'],x['ref'];beta,se,pv,f=map(float,(x['beta'],x['sebeta'],x['pval'],x['af_alt']) )
    else:ch,pos,a,b=x['chromosome'],x['base_pair_location'],x['effect_allele'],x['other_allele'];beta,se,pv,f=map(float,(x['beta'],x['standard_error'],x['p_value'],x['effect_allele_frequency']))
    if not snv(a,b) or not valid(beta,se,pv,f):raise ValueError()
   except (ValueError,KeyError):c['gwas_ineligible']+=1;continue
   index[(ch,pos,frozenset((a,b)))].append((x,a,b,beta,se,pv,f));c['gwas_eligible_rows']+=1
  matched=[]
  for k,qrows in qindex.items():
   if len(qrows)!=1:c['qtl_duplicate_key']+=len(qrows);continue
   candidates=index.get(k,[])
   if len(candidates)!=1:c['gwas_no_unique_key']+=1;continue
   x=qrows[0];orig,a,b,beta,se,pv,f=candidates[0];alt=x['alt'];ref=x['ref']
   if trait=='insomnia' and (a!=alt or b!=ref):c['insomnia_refalt_disagree']+=1;continue
   if a!=alt:beta=-beta;f=1-f
   fq=float(x['maf'])
   if abs(min(f,1-f)-fq)>.15:c['frequency_disagree']+=1;continue
   matched.append({'variant':x['variant'],'rsid':x['rsid'],'chr':x['chromosome'],'pos38':x['position'],'ref':ref,'alt':alt,'qtl_beta':x['beta'],'qtl_se':x['se'],'qtl_p':x['pvalue'],'qtl_minor_freq':fq,'gwas_beta_alt':beta,'gwas_se':se,'gwas_p':pv,'gwas_alt_freq':f});c['matched']+=1
  dest=p/(tissue+'_'+trait+'_matched.tsv');fields=['variant','rsid','chr','pos38','ref','alt','qtl_beta','qtl_se','qtl_p','qtl_minor_freq','gwas_beta_alt','gwas_se','gwas_p','gwas_alt_freq']
  with dest.open('w',newline='') as f:w=csv.DictWriter(f,fields,delimiter='\t');w.writeheader();w.writerows(matched)
  outsummary[tissue+'_'+trait]={'counts':dict(c),'qtl_rows':len(q),'gwas_rows':len(gw),'matched':len(matched),'coverage_eligible_rows':len(matched)/c['qtl_eligible_rows'] if c['qtl_eligible_rows'] else None,'file':dest.name,'sha256':hashlib.sha256(dest.read_bytes()).hexdigest()};print(tissue,trait,len(matched),flush=True)
(p/'match_receipt.json').write_text(json.dumps(outsummary,indent=2)+'\n')
