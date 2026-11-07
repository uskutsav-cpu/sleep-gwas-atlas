import collections,csv,hashlib,json,math,subprocess
from pathlib import Path
from pyliftover import LiftOver
p=Path(__file__).resolve().parent;atlas=p.parents[3]
chain=atlas/'placo_diagnostic/cross_source_followup/hg19ToHg38.over.chain.gz'
sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
assert sha(chain)=='5c0598e500ceb5a78c73086929e8ef993aec309bcafb595139b53d440b125a1d'
lo=LiftOver(str(chain));url='https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90565001-GCST90566000/GCST90565837/GCST90565837.tsv';size=713768673
def read(start,n=8192):
 end=min(size-1,start+n-1)
 raw=subprocess.check_output(['curl','-fsSL','--max-time','90','-D','-','--range',f'{start}-{end}',url])
 headers,data=raw.split(b'\r\n\r\n',1)
 assert b' 206 ' in headers.split(b'\r\n',1)[0] and f'bytes {start}-{end}/{size}'.encode() in headers
 assert len(data)==end-start+1
 return data
head=read(0).split(b'\n',1)[0].decode().rstrip('\r').split('\t')
assert head==['chromosome','base_pair_location','effect_allele','other_allele','beta','standard_error','effect_allele_frequency','p_value','n']
def row_at(offset):
 b=read(offset)
 if offset:b=b.split(b'\n',1)[1]
 else:b=b.split(b'\n',1)[1]
 row=b.split(b'\n',1)[0].decode().split('\t');assert len(row)==len(head)
 return int(row[0]),int(row[1])
def locate(target):
 low=0;high=size-1
 while high-low>8192:
  mid=(low+high)//2
  if row_at(mid)<target:low=mid
  else:high=mid
 return low,high
left=locate((5,103447968));right=locate((5,104447969));start=max(0,left[0]-8192);end=min(size,right[1]+8192);assert end-start<5000000
chunk=read(start,end-start);lines=chunk.split(b'\n')[1:-1];idx=collections.defaultdict(list);counts=collections.Counter();prev=None
for line in lines:
 vals=line.decode().rstrip('\r').split('\t');assert len(vals)==len(head);r=dict(zip(head,vals));k=(int(r['chromosome']),int(r['base_pair_location']))
 if prev is not None:assert k>=prev
 prev=k
 if k[0]!=5 or not 103447968<=k[1]<=104447968:continue
 counts['original_window_rows']+=1;a,b=r['effect_allele'],r['other_allele']
 if len(a)!=1 or len(b)!=1 or a not in 'ACGT' or b not in 'ACGT' or a==b or frozenset((a,b)) in (frozenset('AT'),frozenset('CG')):continue
 m=lo.convert_coordinate('chr5',k[1]-1)
 if len(m)!=1 or m[0][0]!='chr5' or m[0][2]!='+':counts['lift_failure']+=1;continue
 idx[(str(m[0][1]+1),frozenset((a,b)))].append(r)
summary={};issues=[]
for trait in ('insomnia','adhd'):
 path=p.parent/(trait+'_matched.tsv');c=collections.Counter()
 with path.open() as f:
  for h in csv.DictReader(f,delimiter='\t'):
   c['matched_rows']+=1;key=(h['pos38'],frozenset((h['effect_allele'],h['other_allele'])));records=idx.get(key,[])
   if len(records)!=1:
    c['missing_or_duplicate_original_key']+=1;issues.append({'trait':trait,'pos38':h['pos38'],'issue':'missing_or_duplicate','n':len(records)});continue
   r=records[0];same=r['effect_allele']==h['effect_allele'];direction=1 if same else -1
   try:
    checks={'beta':math.isclose(direction*float(r['beta']),float(h['physical_beta']),rel_tol=1e-10,abs_tol=1e-12),'se':math.isclose(float(r['standard_error']),float(h['physical_se']),rel_tol=1e-10,abs_tol=1e-12),'p':math.isclose(float(r['p_value']),float(h['physical_p']),rel_tol=1e-10,abs_tol=1e-300),'frequency':math.isclose(float(r['effect_allele_frequency']) if same else 1-float(r['effect_allele_frequency']),float(h['physical_frequency']),rel_tol=1e-10,abs_tol=1e-12)}
   except ValueError:checks={'numeric':False}
   if all(checks.values()):c['concordant']+=1;c['same_effect_allele' if same else 'swapped_effect_allele']+=1
   else:c['discrepant']+=1;issues.append({'trait':trait,'pos38':h['pos38'],'checks':checks})
 summary[trait]={'matched_file_sha256':sha(path),'counts':dict(c)}
out={'url':url,'scope':'partial original source only','range_start':start,'range_bytes':len(chunk),'chain_sha256':sha(chain),'original_counts':dict(counts),'traits':summary,'issues':issues}
(p/'receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k!='issues'},indent=2));print('issues',len(issues))
