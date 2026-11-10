import csv,gzip,hashlib,io,json,subprocess,pathlib
url='https://zenodo.org/records/17069328/files/mv-Sarcopenia.gwas-ssf.tsv.gz?download=1';expect='d9b5c1e7513e209bc57326264c95c356'
class Hashed:
 def __init__(self,r):self.r=r;self.md5=hashlib.md5();self.n=0
 def read(self,n=-1):
  b=self.r.read(n);self.md5.update(b);self.n+=len(b);return b
 def seekable(self):return False
 def readable(self):return True
proc=subprocess.Popen(['curl','-L','--fail','--max-time','900','-sS',url],stdout=subprocess.PIPE)
with proc.stdout as raw:
 h=Hashed(raw);reader=csv.DictReader(io.TextIOWrapper(gzip.GzipFile(fileobj=h),encoding='utf8',newline=''),delimiter='\t');selected=[];n=0
 for row in reader:
  n+=1
  if n>=123456 and (n-123456)%300000==0:selected.append(dict(row,row_number=n))
 digest=h.md5.hexdigest(); size=h.n
 if proc.wait()!=0 or digest!=expect: raise ValueError('incomplete or checksum mismatch')
with open('work/zenodo_genomewide_systematic_rows.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=['row_number']+reader.fieldnames);w.writeheader();w.writerows(selected)
print(json.dumps({'rows':n,'selected':len(selected),'md5':digest,'bytes':size}))
