import csv,gzip,hashlib,json,math,os,subprocess,sys,time
from collections import Counter
from pathlib import Path
p=Path(__file__).resolve().parent;atlas=p.parents[1];base=atlas/'multiscale_embedding/hfrs_extension/heldout_sleep/group_holdout/bounded_model/mvp_feature_transfer/common_support_run';ref=Path('/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/reference/hm3_grch37_variant_map.tsv.gz');prep=json.loads((base/'prepared_receipt.json').read_text());common_file=base/'MVP.sumstats.gz'
def sha(f):
 h=hashlib.sha256()
 with f.open('rb') as g:
  for x in iter(lambda:g.read(1048576),b''):h.update(x)
 return h.hexdigest()
assert sha(common_file)==prep['prepared_sha256']['MVP'];expected_ref=json.loads((atlas/'replication/mvp_finngen_genomewide/stream_receipt.json').read_text())['reference_sha256'];assert sha(ref)==expected_ref
with gzip.open(common_file,'rt',newline='') as f:common={r['SNP'] for r in csv.DictReader(f,delimiter='\t')}
assert len(common)==991707
refmap={};ambiguous=set();seenrs=set()
with gzip.open(ref,'rt',newline='') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  if r['SNP'] not in common:continue
  a,b=r['A1'],r['A2'];pair=frozenset((a,b))
  if len(a)!=1 or len(b)!=1 or not pair <= set('ACGT') or len(pair)!=2 or pair in (frozenset('AT'),frozenset('CG')):continue
  key=(r['CHR'],r['BP'],''.join(sorted((a,b))))
  if key in refmap:ambiguous.add(key)
  else:refmap[key]=r['SNP']
  assert r['SNP'] not in seenrs;seenrs.add(r['SNP'])
for k in ambiguous:refmap.pop(k,None)
url='https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90565001-GCST90566000/GCST90565837/GCST90565837.tsv';h=hashlib.md5();s=hashlib.sha256();counts=Counter();buf=b'';head=None;records={};duplicate=set();ns=[]
def process(line):
 global head
 cols=line.decode().rstrip('\r').split('\t')
 if head is None:
  head=cols;assert head==['chromosome','base_pair_location','effect_allele','other_allele','beta','standard_error','effect_allele_frequency','p_value','n'];return
 counts['source_rows']+=1
 if len(cols)!=len(head):raise ValueError('malformed row width')
 chr,pos,a,b,beta,se,freq,pval,n=cols
 if len(a)!=1 or len(b)!=1 or a not in 'ACGT' or b not in 'ACGT' or a==b:return
 key=(chr,pos,''.join(sorted((a,b))));rs=refmap.get(key)
 if rs is None:return
 counts['coordinate_allele_rows']+=1
 try:
  beta,se,freq,pval,n=map(float,(beta,se,freq,pval,n))
  if not all(map(math.isfinite,(beta,se,freq,pval,n))) or se<=0 or n<=0 or not 0<freq<1 or not 0<pval<=1:raise ValueError()
  z=beta/se
  if not math.isfinite(z):raise ValueError()
 except (ValueError,OverflowError):counts['invalid_matched']+=1;return
 if rs in records:duplicate.add(rs);counts['duplicate_eligible_rows']+=1;return
 records[rs]=f'{rs}\t{a}\t{b}\t{n:.6f}\t{z:.12g}\n';ns.append(n)
# Resume from the exact byte already hashed after a connection failure. Keep the
# unfinished line in memory; the publisher checksum gates all downstream output.
for attempt in range(11):
 offset=counts['source_bytes']
 cmd=['curl','-fsSL','--connect-timeout','30','--speed-time','120','--speed-limit','1']
 if offset:cmd += ['--range',str(offset)+'-']
 proc=subprocess.Popen(cmd+[url],stdout=subprocess.PIPE)
 while True:
  chunk=proc.stdout.read(1048576)
  if not chunk:break
  h.update(chunk);s.update(chunk);counts['source_bytes']+=len(chunk);buf+=chunk;lines=buf.split(b'\n');buf=lines.pop()
  for line in lines:
   if line:process(line)
 status=proc.wait()
 if counts['source_bytes']>713768673:raise RuntimeError('source exceeded expected size; range may have been ignored')
 if status==0 and counts['source_bytes']==713768673:break
 counts['connection_retries']+=1
 print(f'Connection failed at byte {counts["source_bytes"]}; resuming (attempt {attempt+1}/10).',file=sys.stderr,flush=True)
 if attempt==10:raise RuntimeError('source connection retry limit exceeded')
 time.sleep(min(5*(attempt+1),30))
if buf:process(buf)
assert status==0 and h.hexdigest()=='c10d4dca6efaa4dba2e39ae29add61ae', 'source incomplete or publisher checksum mismatch'
for rs in duplicate:records.pop(rs,None)
counts['duplicate_ids_excluded']=len(duplicate);counts['retained']=len(records);assert len(records)>200000
print('SNP\tA1\tA2\tN\tZ',flush=True)
for line in records.values():print(line,end='')
receipt={'source_url':url,'source_md5':h.hexdigest(),'source_sha256':s.hexdigest(),'common_input_sha256':sha(common_file),'reference_sha256':sha(ref),'reference_ambiguous_keys_excluded':len(ambiguous),'counts':dict(counts),'n_min_including_later_duplicate_exclusions':min(ns),'n_max_including_later_duplicate_exclusions':max(ns),'direction':'source beta / standard_error on reported effect allele','match':'GRCh37 chromosome, position and exact unordered nonpalindromic allele pair; eligible duplicate rsIDs excluded'}
Path(os.environ['STREAM_RECEIPT']).write_text(json.dumps(receipt,indent=2)+'\n')
