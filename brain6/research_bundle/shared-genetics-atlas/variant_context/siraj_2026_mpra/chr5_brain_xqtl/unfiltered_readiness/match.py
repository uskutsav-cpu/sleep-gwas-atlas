import csv,json,math,collections,hashlib
from pathlib import Path
from pyliftover import LiftOver
p=Path(__file__).resolve().parent;atlas=Path('outputs/shared-genetics-atlas');base=Path('/Volumes/Extreme SSD/brain6-work/exploratory-finemap-coloc-v2/inputs');dest=Path('/Volumes/Brain6Work/xqtl-chr5/methylation_candidate');dest.mkdir(exist_ok=True)
chain=atlas/'placo_diagnostic/cross_source_followup/hg19ToHg38.over.chain.gz';lo=LiftOver(str(chain));prov=json.loads(Path('outputs/brain6-alphagenome-pilot/results/provenance.json').read_text())['input_sha256']
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def snv(a,b):return len(a)==len(b)==1 and a in 'ACGT' and b in 'ACGT' and a!=b and frozenset((a,b)) not in (frozenset('AT'),frozenset('CG'))
q=collections.defaultdict(list);qc=collections.Counter()
with (p/'selected_cpg.tsv').open() as f:
 for x in csv.reader(f,delimiter='\t'):
  try:
   ref,alt=x[6:8];beta,se,pv,af=map(float,(x[12],x[13],x[4],x[14]));pos=int(x[2]);
   if not snv(ref,alt) or not all(map(math.isfinite,(beta,se,pv,af))) or se<=0 or not 0<af<1:raise ValueError()
   q[(pos,frozenset((ref,alt)))].append(x);qc['eligible']+=1
  except (ValueError,IndexError):qc['excluded']+=1
summ={};fields=['variant','pos38','ref','alt','qtl_beta','qtl_se','qtl_p','qtl_alt_freq','gwas_beta_alt','gwas_se','gwas_rsid']
for trait in ('insomnia','adhd'):
 name='placo_insomnia__adhd_chr5_103447968_104447968__'+trait+'.tsv';path=base/name;assert sha(path)==prov[name];idx=collections.defaultdict(list);gcounts=collections.Counter()
 with path.open() as f:
  for x in csv.DictReader(f,delimiter='\t'):
   try:
    a,b=x['A1'],x['A2'];beta,se=float(x['BETA']),float(x['SE'])
    if not snv(a,b) or not math.isfinite(beta) or not math.isfinite(se) or se<=0:raise ValueError()
    m=lo.convert_coordinate('chr'+x['CHR'],int(x['POS'])-1)
    if len(m)!=1 or m[0][0]!='chr5' or m[0][2]!='+':gcounts['liftover_fail']+=1;continue
    idx[(m[0][1]+1,frozenset((a,b)))].append(x);gcounts['eligible']+=1
   except (ValueError,KeyError):gcounts['invalid']+=1
 out=[];counts=collections.Counter()
 for key,qs in q.items():
  if len(qs)!=1:counts['qtl_duplicate']+=len(qs);continue
  gs=idx.get(key,[])
  if len(gs)!=1:counts['no_unique_gwas']+=1;continue
  x=qs[0];g=gs[0];alt=x[7];out.append(dict(variant=x[3],pos38=x[2],ref=x[6],alt=alt,qtl_beta=x[12],qtl_se=x[13],qtl_p=x[4],qtl_alt_freq=x[14],gwas_beta_alt=float(g['BETA'])*(1 if g['A1']==alt else -1),gwas_se=g['SE'],gwas_rsid=g['SNP']))
 result=dest/(trait+'.tsv')
 with result.open('w',newline='') as f:
  w=csv.DictWriter(f,fields,delimiter='\t');w.writeheader();w.writerows(out)
 summ[trait]={'gwas_input_sha256':sha(path),'gwas_counts':dict(gcounts),'matching_counts':dict(counts),'matched':len(out),'output_path':str(result),'output_sha256':sha(result)}
(p/'matching_receipt.json').write_text(json.dumps({'qtl_eligibility':dict(qc),'chain_sha256':sha(chain),'pairs':summ,'warning':'selected region; p-value beta/SE are based on t370, Gaussian ABFs are approximate; GWAS allele frequency not available'},indent=2)+'\n');print(json.dumps(summ,indent=2))
