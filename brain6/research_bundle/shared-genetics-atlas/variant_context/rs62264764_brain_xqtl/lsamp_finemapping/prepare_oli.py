import collections,csv,hashlib,json,math
from pathlib import Path
from pyliftover import LiftOver
import pysam
p=Path(__file__).resolve().parent;atlas=p.parents[2];out=Path('/Volumes/Brain6Work/xqtl-chr3/lsamp_matched');out.mkdir(exist_ok=True)
base=Path('/Volumes/Brain6Work/xqtl-chr3');ext=Path('/Volumes/Extreme SSD/brain6-work/exploratory-finemap-coloc-v2/inputs');prov=json.loads((p.parents[3]/'brain6-alphagenome-pilot/results/provenance.json').read_text())['input_sha256'];chain=atlas/'placo_diagnostic/cross_source_followup/hg19ToHg38.over.chain.gz';lo=LiftOver(str(chain)); lower,upper=117423158,118423158
expected={'oli.bed.gz':'363b2c8ed2575bf7d4da68c5d9bc8762','oli.bed.gz.tbi':'d1278d1ca182fc08b779b18fa224321e'}
def digest(path,method='sha256'):
 h=hashlib.new(method)
 with path.open('rb') as f:
  for block in iter(lambda:f.read(4*1024*1024),b''):h.update(block)
 return h.hexdigest()
for name,md5 in expected.items():assert digest(base/name,'md5')==md5
for b,pos in [(117142005,lower),(118142005,upper)]:
 m=lo.convert_coordinate('chr3',b-1);assert len(m)==1 and m[0][0]=='chr3' and m[0][2]=='+' and m[0][1]+1==pos
snv=lambda a,b:len(a)==len(b)==1 and a in 'ACGT' and b in 'ACGT' and a!=b and frozenset((a,b)) not in (frozenset('AT'),frozenset('CG'))
gwas={};source={}
for t in ('insomnia','adhd'):
 name='placo_insomnia__adhd_chr3_117142005_118142005__'+t+'.tsv'; path=ext/name;assert digest(path)==prov[name]
 ix=collections.defaultdict(list);ct=collections.Counter()
 with path.open() as f:
  for r in csv.DictReader(f,delimiter='\t'):
   ct['total']+=1
   try:
    a,b=r['A1'],r['A2'];beta,se=map(float,(r['BETA'],r['SE']))
    if not snv(a,b) or not math.isfinite(beta) or not math.isfinite(se) or se<=0:raise ValueError()
    m=lo.convert_coordinate('chr'+r['CHR'],int(r['POS'])-1)
    if len(m)!=1 or m[0][0]!='chr3' or m[0][2]!='+' or not lower<=m[0][1]+1<=upper:ct['lift_or_window_fail']+=1;continue
    ix[(m[0][1]+1,frozenset((a,b)))].append(r);ct['eligible']+=1
   except (ValueError,KeyError):ct['invalid']+=1
 gwas[t]=ix;source[t]={'path':str(path),'sha256':digest(path),'counts':dict(ct)}
summ={};cols=['variant','pos38','ref','alt','qtl_beta','qtl_se','qtl_p','qtl_alt_freq','gwas_rsid','gwas_beta_alt','gwas_se']
for ctx,gene in [('oli','ENSG00000185565')]:
 tab=pysam.TabixFile(str(base/(ctx+'.bed.gz')));header=list(tab.header)[0].lstrip('#').split('\t'); q=collections.defaultdict(list);ct=collections.Counter();exact_variant=None
 for line in tab.fetch('chr3',lower-1,upper):
  r=dict(zip(header,line.split('\t')))
  if r['target_ensembl_id']!=gene:continue
  ct['target_rows']+=1
  try:
   ref,alt=r['ref'],r['alt'];b,se,pv,af=map(float,(r['beta_non_ref'],r['beta_se_non_ref'],r['pval'],r['non_ref_af']));pos=int(r['chromEnd'])
   if not snv(ref,alt) or not all(map(math.isfinite,(b,se,pv,af))) or se<=0 or not(0<=pv<=1 and 0<af<1):raise ValueError()
   if r['variant_id']=='chr3:117920728:G:A':
    exact_variant={'beta':b,'se':se,'p':pv}
   q[(pos,frozenset((ref,alt)))].append(r);ct['eligible']+=1
  except (ValueError,KeyError):ct['invalid']+=1
 tab.close(); assert exact_variant is not None
 for t in ('insomnia','adhd'):
  count=collections.Counter();rows=[]
  for key,qr in q.items():
   if len(qr)!=1:count['duplicate_qtl_key']+=len(qr);continue
   gr=gwas[t].get(key,[])
   if len(gr)!=1:count['gwas_no_unique_key']+=1;continue
   x=qr[0];g=gr[0];beta=float(g['BETA'])*(1 if g['A1']==x['alt'] else -1)
   rows.append({'variant':x['variant_id'],'pos38':key[0],'ref':x['ref'],'alt':x['alt'],'qtl_beta':x['beta_non_ref'],'qtl_se':x['beta_se_non_ref'],'qtl_p':x['pval'],'qtl_alt_freq':x['non_ref_af'],'gwas_rsid':g['SNP'],'gwas_beta_alt':beta,'gwas_se':g['SE']});count['matched']+=1
  dest=out/(ctx+'_'+t+'.tsv')
  with dest.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=cols,delimiter='\t');w.writeheader();w.writerows(rows)
  summ[ctx+'_'+t]={'gene':gene,'qtl_counts':dict(ct),'match_counts':dict(count),'file':str(dest),'sha256':digest(dest),'exact_rs62264764_nominal':exact_variant}
  print(ctx,t,len(rows),flush=True)
receipt={'source':source,'chain_sha256':digest(chain),'qtl_source':{'dataset':'NG00184.v1','file_md5':expected,'file_sha256':{name:digest(base/name) for name in expected},'context':'ROSMAP_CUIMC1_2_MIT','window_grch38':[lower,upper]},'pairs':summ,'note':'No GWAS source frequencies; selected region and incomplete cis/TAD support; matches do not prove phenotype causality.'}
(p/'matching_receipt_oli.json').write_text(json.dumps(receipt,indent=2)+'\n')
