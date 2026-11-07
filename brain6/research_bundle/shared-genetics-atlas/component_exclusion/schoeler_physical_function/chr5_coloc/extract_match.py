import collections,csv,gzip,hashlib,json,math,subprocess,zlib
from pathlib import Path
from pyliftover import LiftOver
import pysam
p=Path(__file__).resolve().parent;atlas=p.parents[2];workspace=atlas.parents[1]
source=Path('/Volumes/Extreme SSD/brain6-work/exploratory-finemap-coloc-v2/inputs')
provenance=atlas.parent/'brain6-alphagenome-pilot/results/provenance.json'
chain=atlas/'placo_diagnostic/cross_source_followup/hg19ToHg38.over.chain.gz'
sha=lambda q:hashlib.sha256(q.read_bytes()).hexdigest()
assert sha(chain)=='5c0598e500ceb5a78c73086929e8ef993aec309bcafb595139b53d440b125a1d'
prov=json.loads(provenance.read_text())['input_sha256'];lo=LiftOver(str(chain))
url='https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90565001-GCST90566000/GCST90565837/harmonised/GCST90565837.h.tsv.gz'
chunk=subprocess.check_output(['curl','-fsSL','--range','0-65535','--max-time','45',url]);head=zlib.decompressobj(31).decompress(chunk).split(b'\n',1)[0].decode().rstrip('\r').split('\t')
assert head==['chromosome','base_pair_location','effect_allele','other_allele','beta','standard_error','effect_allele_frequency','p_value','rsid','n','hm_coordinate_conversion','hm_code','variant_id']
idxpath=p.parent/'brain6_loci/GCST90565837.h.tsv.gz.tbi';f=pysam.TabixFile(url,index=str(idxpath));physical=collections.defaultdict(list);count=collections.Counter()
def valid(a,b):return len(a)==len(b)==1 and a in 'ACGT' and b in 'ACGT' and a!=b and frozenset((a,b)) not in (frozenset('AT'),frozenset('CG'))
for line in f.fetch('5',104112266,105112267):
 r=dict(zip(head,line.split('\t')));count['physical_rows']+=1;a,b=r['effect_allele'],r['other_allele']
 try:
  beta,se,freq,pv,n=map(float,(r['beta'],r['standard_error'],r['effect_allele_frequency'],r['p_value'],r['n']))
  if not valid(a,b) or not all(map(math.isfinite,(beta,se,freq,pv,n))) or se<=0 or not 0<freq<1 or not 0<pv<=1 or n<=0:raise ValueError()
  physical[(r['base_pair_location'],frozenset((a,b)))].append(r);count['physical_eligible']+=1
 except (ValueError,KeyError):count['physical_invalid']+=1
results={}
for trait in ('insomnia','adhd'):
 name='placo_insomnia__adhd_chr5_103447968_104447968__'+trait+'.tsv';file=source/name;h=sha(file);assert h==prov[name]
 src=collections.defaultdict(list);counts=collections.Counter()
 with file.open() as fh:
  for r in csv.DictReader(fh,delimiter='\t'):
   counts['source_rows']+=1
   try:
    a,b=r['A1'],r['A2'];beta,se=float(r['BETA']),float(r['SE']);m=lo.convert_coordinate('chr'+r['CHR'],int(r['POS'])-1)
    if not valid(a,b) or not math.isfinite(beta) or not math.isfinite(se) or se<=0 or len(m)!=1 or m[0][0]!='chr5' or m[0][2]!='+':raise ValueError()
    pos=m[0][1]+1
    if not 104112267<=pos<=105112267:raise ValueError()
    src[(str(pos),frozenset((a,b)))].append(r);counts['source_eligible']+=1
   except (ValueError,KeyError):counts['source_invalid']+=1
 out=[]
 for key,rs in src.items():
  q=physical.get(key,[])
  if len(rs)!=1 or len(q)!=1:counts['no_unique_match']+=1;continue
  r=rs[0];t=q[0];sign=1 if r['A1']==t['effect_allele'] else -1
  out.append({'rsid_source':r['SNP'],'pos38':key[0],'effect_allele':t['effect_allele'],'other_allele':t['other_allele'],'trait_beta':sign*float(r['BETA']),'trait_se':float(r['SE']),'physical_beta':float(t['beta']),'physical_se':float(t['standard_error']),'physical_p':float(t['p_value']),'physical_frequency':float(t['effect_allele_frequency']),'physical_rsid':t['rsid']})
 counts['matched']=len(out);out.sort(key=lambda x:int(x['pos38']))
 path=p/(trait+'_matched.tsv')
 if out:
  with path.open('w',newline='') as fh:w=csv.DictWriter(fh,fieldnames=list(out[0]),delimiter='\t');w.writeheader();w.writerows(out)
 results[trait]={'source_sha256':h,'counts':dict(counts),'matched_sha256':sha(path) if out else None}
receipt={'url':url,'partial_indexed_query':True,'index_sha256':sha(idxpath),'chain_sha256':sha(chain),'region_grch38_inclusive':[104112267,105112267],'physical_counts':dict(count),'results':results}
(p/'match_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
