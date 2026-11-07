import pathlib,csv,json,collections,hashlib,math
from pyliftover import LiftOver
p=pathlib.Path(__file__).resolve().parent
external=pathlib.Path('/Volumes/Extreme SSD/brain6-work/exploratory-finemap-coloc-v2/inputs')
provenance=p.parents[2]/'brain6-alphagenome-pilot/results/provenance.json'
chain=p.parents[1]/'placo_diagnostic/cross_source_followup/hg19ToHg38.over.chain.gz'
prov=json.loads(provenance.read_text())['input_sha256'];lo=LiftOver(str(chain))
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def read(path):
 with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def snv(a,b):return len(a)==len(b)==1 and a in 'ACGT' and b in 'ACGT' and a!=b and frozenset((a,b)) not in (frozenset('AT'),frozenset('CG'))
source={};indices={};gcounts={}
for trait in ('insomnia','adhd'):
 name='placo_insomnia__adhd_chr3_117142005_118142005__'+trait+'.tsv';path=external/name;h=sha(path);assert prov[name]==h
 data=read(path);idx=collections.defaultdict(list);c=collections.Counter()
 for x in data:
  try:
   beta,se=float(x['BETA']),float(x['SE']);a,b=x['A1'],x['A2']
   if not snv(a,b) or not math.isfinite(beta) or not math.isfinite(se) or se<=0:raise ValueError()
   m=lo.convert_coordinate('chr'+x['CHR'],int(x['POS'])-1)
   if len(m)!=1 or m[0][0]!='chr3' or m[0][2]!='+':c['lift_failed']+=1;continue
   pos=m[0][1]+1
   if not 117423158<=pos<=118423158:c['outside_window']+=1;continue
   idx[(str(pos),frozenset((a,b)))].append(x);c['eligible_rows']+=1
  except (ValueError,KeyError):c['ineligible_rows']+=1
 indices[trait]=idx;source[trait]={'input_path':str(path),'sha256':h,'total_rows':len(data),'counts':dict(c)}
summary={};fields=['variant','qtl_rsid','gwas_rsid','pos38','ref','alt','qtl_beta','qtl_se','qtl_p','qtl_maf','gwas_beta_alt','gwas_se','gwas_p_from_beta_se']
for tissue in ('Adrenal_Gland','Brain_Frontal_Cortex_BA9'):
 for gene in ('ENSG00000185565','ENSG00000239268'):
  path=p/(tissue+'_'+gene+'.tsv');rows=read(path);qidx=collections.defaultdict(list);qc=collections.Counter()
  for x in rows:
   try:
    a,b=x['ref'],x['alt'];beta,se,pv,maf=map(float,(x['beta'],x['se'],x['pvalue'],x['maf']))
    if not snv(a,b) or not(all(map(math.isfinite,(beta,se,pv,maf))) and se>0 and 0<=pv<=1 and 0<maf<=.5):raise ValueError()
   except ValueError:qc['ineligible_rows']+=1;continue
   qc['eligible_rows']+=1;qidx[(x['position'],frozenset((a,b)))].append(x)
  for trait in ('insomnia','adhd'):
   c=qc.copy();matched=[]
   for key,records in qidx.items():
    if len(records)!=1:c['duplicate_qtl_key']+=len(records);continue
    src=indices[trait].get(key,[])
    if len(src)!=1:c['gwas_no_unique_key']+=1;continue
    q=records[0];g=src[0];beta=float(g['BETA'])*(1 if g['A1']==q['alt'] else -1);se=float(g['SE']);pv=math.erfc(abs(beta/se)/math.sqrt(2))
    matched.append({'variant':q['variant'],'qtl_rsid':q['rsid'],'gwas_rsid':g['SNP'],'pos38':q['position'],'ref':q['ref'],'alt':q['alt'],'qtl_beta':q['beta'],'qtl_se':q['se'],'qtl_p':q['pvalue'],'qtl_maf':q['maf'],'gwas_beta_alt':beta,'gwas_se':se,'gwas_p_from_beta_se':pv});c['matched']+=1
   dest=p/(tissue+'_'+gene+'_'+trait+'_matched.tsv')
   with dest.open('w',newline='') as f:w=csv.DictWriter(f,fields,delimiter='\t');w.writeheader();w.writerows(matched)
   summary[tissue+'_'+gene+'_'+trait]={'counts':dict(c),'qtl_total':len(rows),'matched':len(matched),'eligible_coverage':len(matched)/qc['eligible_rows'] if qc['eligible_rows'] else None,'sha256':sha(dest),'file':dest.name};print(tissue,gene,trait,len(matched),flush=True)
(p/'match_receipt.json').write_text(json.dumps({'source':source,'chain_sha256':sha(chain),'pairs':summary},indent=2)+'\n')
