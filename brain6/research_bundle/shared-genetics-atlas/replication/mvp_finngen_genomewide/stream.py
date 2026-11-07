import csv,gzip,hashlib,json,math,os,sys
from collections import Counter
from pathlib import Path
sys.path.insert(0,'/Users/swethasunilkumar/Library/Python/3.13/lib/python/site-packages')
from pyliftover import LiftOver
from scipy.special import ndtri
p=Path(__file__).resolve().parent
ref=Path('/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/reference/hm3_grch37_variant_map.tsv.gz'); chain=Path('work/hg19ToHg38.over.chain.gz'); adhd=Path('/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/Sep1-local-dependencies/data/munged/track_b_finngen_r13_F5_ADHD.sumstats.gz'); mvp=Path('/Volumes/Extreme SSD/sleep-frailty-full-autonomous-v5-20261009/frailty_paper/full_autonomous_research_v5/genetics/new_public_validation/raw/GCST90475826.tsv.gz')
expected={adhd:'24d4d72ebdf1f83e7af4dd5e46aa5208266e121f14bd3d5a26217a87f861f59b',mvp:'46b62349ba0d98b22bf302be973bb54a80707cfb02ce953426f8224910761835',chain:'5c0598e500ceb5a78c73086929e8ef993aec309bcafb595139b53d440b125a1d'}
def sha(f):
 h=hashlib.sha256()
 with f.open('rb') as q:
  while c:=q.read(4*1024*1024):h.update(c)
 return h.hexdigest()
for f,d in expected.items():
 if sha(f)!=d:raise ValueError('input hash mismatch '+str(f))
with gzip.open(adhd,'rt') as f:rows=list(csv.DictReader(f,delimiter='\t'))
common={x['SNP']:frozenset((x['A1'],x['A2'])) for x in rows}; assert len(common)==len(rows)==1156825
lift=LiftOver(str(chain));refmap={};count=Counter(); complement=str.maketrans('ACGT','TGCA')
with gzip.open(ref,'rt') as f:
 for x in csv.DictReader(f,delimiter='\t'):
  if x['SNP'] not in common:continue
  pair=frozenset((x['A1'],x['A2']))
  if pair!=common[x['SNP']] or pair in (frozenset('AT'),frozenset('CG')):count['reference_pair_or_palindrome']+=1;continue
  r=lift.convert_coordinate('chr'+x['CHR'],int(x['BP'])-1)
  if len(r)!=1 or r[0][0]!='chr'+x['CHR'] or r[0][2] not in ('+','-'):count['liftover_excluded']+=1;continue
  if x['SNP'] in refmap:raise ValueError('duplicate reference SNP')
  strand=r[0][2];mapped=pair if strand=='+' else frozenset(v.translate(complement) for v in pair)
  refmap[x['SNP']]=(x['CHR'],int(r[0][1])+1,mapped,strand)
records={}; duplicate_ids=set();Nmin=float('inf');Nmax=0
with gzip.open(mvp,'rt') as f:
 for x in csv.DictReader(f,delimiter='\t'):
  count['rows']+=1;rs=x['rsid']; v=refmap.get(rs)
  if v is None:continue
  count['rsid_matched']+=1;ea=x['effect_allele'];oa=x['other_allele']
  try:
   pos=int(x['base_pair_location']);pair=frozenset((ea,oa))
   if x['chromosome']!=v[0] or pos!=v[1] or pair!=v[2]:count['position_or_allele_mismatch']+=1;continue
   o=float(x['odds_ratio']);pval=float(x['p_value']);nc=float(x['num_cases']);nu=float(x['num_controls']);freq=float(x['effect_allele_frequency'])
   if not all(map(math.isfinite,[o,pval,nc,nu,freq])) or o<=0 or not 0<pval<=1 or nc<=0 or nu<=0 or not 0<freq<1:count['invalid']+=1;continue
   if o==1:count['rounded_OR_no_direction']+=1;continue
   z=(-float(ndtri(pval/2)))*(1 if o>1 else -1);n=4*nc*nu/(nc+nu)
   if not math.isfinite(z) or not math.isfinite(n) or z==0:count['invalid_z_or_n']+=1;continue
   if rs in records:count['duplicate_rows']+=1;duplicate_ids.add(rs);continue
   Nmin=min(Nmin,n);Nmax=max(Nmax,n)
   if v[3]=='-':ea=ea.translate(complement);oa=oa.translate(complement)
   records[rs]=f'{rs}\t{ea}\t{oa}\t{n:.6f}\t{z:.12g}'
  except (ValueError,OverflowError):count['invalid']+=1
count['duplicate_ids_excluded']=len(duplicate_ids);count['retained']=len(records)-len(duplicate_ids)
assert count['rows']==19703815 and count['retained']>200000
print('SNP\tA1\tA2\tN\tZ',flush=True)
for rs,line in records.items():
 if rs not in duplicate_ids:print(line)
Path(os.environ.get('STREAM_RECEIPT', str(p/'stream_receipt.json'))).write_text(json.dumps({'input_sha256':{str(k):v for k,v in expected.items()},'reference_sha256':sha(ref),'reference_lifted':len(refmap),'counts':dict(count),'n_effective_min':Nmin,'n_effective_max':Nmax,'p_to_z':'signed inverse normal from source p; direction from rounded odds ratio'},indent=2)+'\n')
