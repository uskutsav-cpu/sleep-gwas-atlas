import collections,hashlib,json,math,subprocess,zlib
from pathlib import Path
p=Path(__file__).resolve().parent;u='https://zenodo.org/api/records/17069328/files/mv-Sarcopenia.gwas-ssf.tsv.gz/content';proc=subprocess.Popen(['curl','-fsSL','--retry','2',u],stdout=subprocess.PIPE);md=hashlib.md5();decomp=zlib.decompressobj(31);boundary=True;buf=b'';members=0;count=0;bad=collections.Counter();examples={};header=None

def row(line):
 global count,header
 a=line.decode().rstrip('\r').split('\t')
 if header is None:header=a;return
 count+=1;d=dict(zip(header,a));errors=[]
 for key in ('beta','standard_error','p_value','effect_allele_frequency','sample_size'):
  try:
   v=float(d[key]);ok=math.isfinite(v) and (v>0 if key in ('standard_error','sample_size') else 0<=v<=1 if key in ('p_value','effect_allele_frequency') else True)
  except (ValueError,OverflowError):ok=False
  if not ok:errors.append(key)
 if errors:
  bad['any']+=1
  for k in errors:
   bad[k]+=1
   if k not in examples:examples[k]=d[k]
while chunk:=proc.stdout.read(2**20):
 md.update(chunk);remain=chunk
 while remain:
  text=decomp.decompress(remain);remain=decomp.unused_data;boundary=decomp.eof;buf+=text;lines=buf.split(b'\n');buf=lines.pop()
  for line in lines:
   if line:row(line)
  if boundary:members+=1;decomp=zlib.decompressobj(31)
  else:break
if buf:row(buf)
assert proc.wait()==0 and md.hexdigest()=='d9b5c1e7513e209bc57326264c95c356' and boundary and count==7319145 and bad['any']==7396
x={'rows':count,'md5':md.hexdigest(),'gzip_members_crc_checked':members,'numeric_failures':dict(bad),'example_values':examples};(p/'invalid_receipt.json').write_text(json.dumps(x,indent=2)+'\n');print(json.dumps(x,indent=2))
