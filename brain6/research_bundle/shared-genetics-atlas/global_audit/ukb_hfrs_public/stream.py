import collections,csv,gzip,hashlib,json,math,subprocess,zlib
from pathlib import Path

p=Path(__file__).resolve().parent
common=Path('/Volumes/Brain6Work/physical-function-source-sensitivity-20261028/questionnaire.sumstats.gz')
ref=Path('/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/reference/hm3_grch37_variant_map.tsv.gz')
target=Path('/Volumes/Brain6Work/ukb_hfrs_public_20261028.sumstats.gz')
url='https://ndownloader.figshare.com/files/61648825';expected='4eafec0dc6a5e07612d22873ae45bb79'
def sha(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
assert sha(common)=='c66d2fe5eb78edb5613a8b6d06259200e798e2ab065f520c518a1fca1d1c7d50'
assert sha(ref)=='6775a7a0d3ca90dc74e472180b1d77103bc238129c4f969358f46307e5c306b4'
want={}
with gzip.open(common,'rt') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  assert r['SNP'] not in want;want[r['SNP']]=(r['A1'],r['A2'])
refs={}
with gzip.open(ref,'rt') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  k=r['SNP']
  if k in want:
   assert k not in refs;refs[k]=(r['CHR'],r['BP'],r['A1'],r['A2'])
assert len(refs)==len(want)
process=subprocess.Popen(['curl','-fsSL','--retry','2',url],stdout=subprocess.PIPE)
md5=hashlib.md5();digest=hashlib.sha256();z=zlib.decompressobj(31);buffer=b'';size=0;header=None
counts=collections.Counter();n_counts=collections.Counter();candidates={};duplicates=set();comp=str.maketrans('ACGT','TGCA')
def consume(line):
 global header
 if header is None:
  header=line.decode().split('\t');assert header==['CHR','SNP','POS','A1','A2','N','AF1','BETA','SE','P','INFO'];return
 v=line.decode().split('\t');assert len(v)==11
 ch,k,pos,a,b,n,af,beta,se,pval,info=v;counts['source_rows']+=1;n_counts[n]+=1
 if k not in want:return
 counts['common_id_rows']+=1
 rc,rp,ra,rb=refs[k]
 if (ch,pos)!=(rc,rp):counts['coordinate_mismatch']+=1;return
 counts['coordinate_match_rows']+=1
 if len(a)!=1 or len(b)!=1 or a not in 'ACGT' or b not in 'ACGT' or a==b or {a,b} in ({'A','T'},{'C','G'}):counts['non_snv_or_palindromic']+=1;return
 if {a,b}!={ra,rb} and {a.translate(comp),b.translate(comp)}!={ra,rb}:counts['reference_allele_mismatch']+=1;return
 if {a,b}!=set(want[k]) and {a.translate(comp),b.translate(comp)}!=set(want[k]):counts['common_allele_mismatch']+=1;return
 counts['matched_allele_rows']+=1
 if k in candidates or k in duplicates:duplicates.add(k);candidates.pop(k,None);return
 try:
  nn,aa,be,ss,pp,ii=map(float,(n,af,beta,se,pval,info))
  if not all(map(math.isfinite,(nn,aa,be,ss,pp,ii))) or nn<=0 or ss<=0 or not(0<=aa<=1 and 0<=pp<=1 and 0<=ii<=1):raise ValueError()
 except ValueError:counts['invalid_numeric']+=1;candidates[k]=None;return
 if min(aa,1-aa)<0.01 or ii<0.9:counts['failed_maf_or_info']+=1;candidates[k]=None;return
 zv=be/ss;approx=math.erfc(abs(zv)/math.sqrt(2))
 if pp>0 and approx>0:
  counts['p_compared']+=1
  if abs(math.log10(pp)-math.log10(approx))>.02:counts['p_log10_difference_gt_0.02']+=1
 candidates[k]=(k,a,b,n,format(zv,'.12g'))
while chunk:=process.stdout.read(1048576):
 size+=len(chunk);md5.update(chunk);digest.update(chunk);buffer+=z.decompress(chunk)
 lines=buffer.split(b'\n');buffer=lines.pop()
 for line in lines:
  if line:consume(line.rstrip(b'\r'))
if buffer:consume(buffer.rstrip(b'\r'))
code=process.wait();assert code==0 and z.eof and size==280558723 and md5.hexdigest()==expected,(code,z.eof,size,md5.hexdigest())
counts['duplicate_matching_ids_excluded']=len(duplicates)
selected={k:v for k,v in candidates.items() if v is not None and k not in duplicates}
counts['selected']=len(selected)
with gzip.open(target,'xt') as f:
 w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['SNP','A1','A2','N','Z'])
 for k in want:
  if k in selected:w.writerow(selected[k])
out={'source_url':url,'source_name':'UKB_hfrs_gwas.txt.gz','source_size':size,'source_md5':md5.hexdigest(),
     'publisher_md5_matches':True,'source_sha256':digest.hexdigest(),'gzip_complete':True,'n_counts':dict(n_counts),
     'published_n_description':407463,'counts':dict(counts),'common_source_sha256':sha(common),
     'reference_sha256':sha(ref),'output_path':str(target),'output_sha256':sha(target)}
(p/'stream_receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
