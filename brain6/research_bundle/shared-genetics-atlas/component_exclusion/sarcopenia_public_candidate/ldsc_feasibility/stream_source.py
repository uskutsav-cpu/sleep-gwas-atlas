import csv,gzip,hashlib,json,math,subprocess,sys,zlib,os
from pathlib import Path
p=Path(__file__).resolve().parent
prior=p.parents[1]/'physical_function_ldsc/insomnia_source_sensitivity/prepare_receipt.json';pr=json.loads(prior.read_text());common={}
with gzip.open(pr['output_paths']['questionnaire'],'rt') as f:
 for r in csv.DictReader(f,delimiter='\t'):common[r['SNP']]=frozenset((r['A1'],r['A2']))
assert len(common)==1020637
ref={};path='/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/reference/hm3_grch37_variant_map.tsv.gz'
with gzip.open(path,'rt') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  if r['SNP'] in common:
   assert r['SNP'] not in ref;ref[r['SNP']]=(r['CHR'],r['BP'],frozenset((r['A1'],r['A2'])))
u='https://zenodo.org/api/records/17069328/files/mv-Sarcopenia.gwas-ssf.tsv.gz/content';c=subprocess.Popen(['curl','-fsSL','--retry','2',u],stdout=subprocess.PIPE);h=hashlib.md5();dec=zlib.decompressobj(31);buf=b'';done=True;members=0;size=0;rows=0;kept=0;seen=set();duplicates=0;bad=0;head=None
print('SNP\tA1\tA2\tN\tZ',flush=True)
def process(line):
 global head,rows,kept,duplicates,bad
 cols=line.decode().rstrip('\r').split('\t')
 if head is None:head=cols;return
 rows+=1
 ch,pos,a,b,beta,se,pv,freq,rs,n=cols
 v=ref.get(rs)
 if v is None or v[0]!=ch or v[1]!=pos or v[2]!=frozenset((a,b)) or v[2]!=common[rs] or v[2] in (frozenset('AT'),frozenset('CG')):return
 try:
  beta,se,pv,freq,n=map(float,(beta,se,pv,freq,n));assert all(map(math.isfinite,(beta,se,pv,freq,n))) and se>0 and 0<=pv<=1 and 0<freq<1 and n==522632
 except (ValueError,AssertionError,OverflowError):bad+=1;return
 if rs in seen:duplicates+=1;return
 seen.add(rs);kept+=1;print(f'{rs}\t{a}\t{b}\t{n:g}\t{beta/se:.12g}')
while chunk:=c.stdout.read(2**20):
 h.update(chunk);size+=len(chunk);rem=chunk
 while rem:
  text=dec.decompress(rem);rem=dec.unused_data;done=dec.eof;buf+=text;lines=buf.split(b'\n');buf=lines.pop()
  for line in lines:
   if line:process(line)
  if done:members+=1;dec=zlib.decompressobj(31)
  else:break
if buf:process(buf)
assert c.wait()==0 and done and size==269614908 and h.hexdigest()=='d9b5c1e7513e209bc57326264c95c356'
out={'source_bytes':size,'source_md5':h.hexdigest(),'rows':rows,'retained':kept,'duplicate_eligible_rows':duplicates,'excluded_invalid_common_rows':bad,'gzip_members':members,'n_assumed':522632,'warning':'Factor inputs, N meaning, factor sign and overlap unresolved.'}
(p/('stream_receipt_'+os.environ['STREAM_TAG']+'.json')).write_text(json.dumps(out,indent=2)+'\n')
