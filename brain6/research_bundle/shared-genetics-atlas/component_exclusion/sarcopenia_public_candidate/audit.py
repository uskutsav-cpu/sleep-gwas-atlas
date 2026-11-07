import csv,gzip,hashlib,json,math,subprocess,zlib
from pathlib import Path
p=Path(__file__).resolve().parent;atlas=p.parents[1]
ref=Path('/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/reference/hm3_grch37_variant_map.tsv.gz')
refs={}
with gzip.open(ref,'rt') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  if r['SNP'] in refs:refs[r['SNP']]=None
  else:refs[r['SNP']]=(r['CHR'],r['BP'],frozenset((r['A1'],r['A2'])))
assert len(refs)>1000000
sleep=[]
with (atlas/'replication/six_lead_joint/six_leads.tsv').open() as f:
 sleep=list(csv.DictReader(f,delimiter='\t'))
assert len(sleep)==6;targets={r['variant'] for r in sleep}|{'rs62264764'}
u='https://zenodo.org/api/records/17069328/files/mv-Sarcopenia.gwas-ssf.tsv.gz/content';proc=subprocess.Popen(['curl','-fsSL','--retry','2',u],stdout=subprocess.PIPE)
md=hashlib.md5();sha=hashlib.sha256();decoder=zlib.decompressobj(31);boundary=True;buffer=b'';size=rows=members=0;header=None;matches=[];stats={'reference_id_rows':0,'grch37_pos_match':0,'grch37_pos_mismatch':0,'ref_pair_match':0,'ref_pair_mismatch':0,'invalid_numeric':0,'target_rows':0};Ns=set();q=[]
def row(line):
 global header,rows
 s=line.decode().rstrip('\r').split('\t')
 if header is None:
  header=s;assert header==['chromosome','base_pair_location','effect_allele','other_allele','beta','standard_error','p_value','effect_allele_frequency','rsid','sample_size'];return
 assert len(s)==len(header);rows+=1;ch,pos,a,b,beta,se,pv,freq,rs,n=s
 try:
  beta,se,pv,freq,n=map(float,(beta,se,pv,freq,n));Ns.add(n)
  assert all(map(math.isfinite,(beta,se,pv,freq,n))) and se>0 and 0<=pv<=1 and 0<=freq<=1 and n>0
  if len(q)<10000 and pv>0:q.append(abs(math.erfc(abs(beta/se)/math.sqrt(2))-pv))
 except (ValueError,AssertionError,OverflowError):stats['invalid_numeric']+=1
 if rs in targets:matches.append(dict(zip(header,s)));stats['target_rows']+=1
 v=refs.get(rs)
 if v is not None:
  stats['reference_id_rows']+=1
  if (ch,pos)==v[:2]:
   stats['grch37_pos_match']+=1
   if frozenset((a,b))==v[2]:stats['ref_pair_match']+=1
   else:stats['ref_pair_mismatch']+=1
  else:stats['grch37_pos_mismatch']+=1
while chunk:=proc.stdout.read(2**20):
 size+=len(chunk);md.update(chunk);sha.update(chunk);remaining=chunk
 while remaining:
  data=decoder.decompress(remaining);remaining=decoder.unused_data;boundary=decoder.eof;buffer+=data;lines=buffer.split(b'\n');buffer=lines.pop()
  for line in lines:
   if line:row(line)
  if decoder.eof:members+=1;decoder=zlib.decompressobj(31)
  else:break
 if size%(64*2**20)==0:print('streamed',size,'rows',rows,flush=True)
if buffer:row(buffer)
assert proc.wait()==0 and size==269614908 and md.hexdigest()=='d9b5c1e7513e209bc57326264c95c356' and boundary
(p/'selected_source_rows.json').write_text(json.dumps(matches,indent=2)+'\n')
receipt={'url':u,'bytes':size,'md5':md.hexdigest(),'sha256':sha.hexdigest(),'gzip_members_crc_checked':members,'source_rows':rows,'counts':stats,'distinct_N':sorted(Ns),'max_abs_p_discrepancy_first_10000':max(q) if q else None,'reference':str(ref),'selected_ids':sorted(targets),'source_model':'depositor describes genomic SEM; full input and N provenance unresolved'}
(p/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
