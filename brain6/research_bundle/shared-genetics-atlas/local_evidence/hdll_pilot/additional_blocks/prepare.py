import csv, gzip, hashlib, io, json, shutil, zipfile
from pathlib import Path
import requests

out=Path(__file__).resolve().parent
work=Path('/Volumes/Brain6Work/hdll-additional-20261028')
work.mkdir(exist_ok=False)
base='https://zenodo.org/api/records/14825987/files/'
b=requests.get(base+'bimfile.zip/content',timeout=90);b.raise_for_status();data=b.content
assert len(data)==9673272 and hashlib.md5(data).hexdigest()=='608f037e4fc4d69f0ee9231c7fc005bb'
leads={'rs62264767':3,'rs2431108':5,'rs1476535':7};selected={}
with zipfile.ZipFile(io.BytesIO(data)) as z:
 for lead,chrom in leads.items():
  hits=[]
  for item in z.infolist():
   if not item.filename.startswith(f'bimfile/ukb_chr{chrom}.') or not item.filename.endswith('.bim'): continue
   bb=z.read(item);rows=[l.split() for l in bb.decode().splitlines()]
   if any(r[1]==lead for r in rows): hits.append((item,bb,rows))
  if len(hits)!=1:
   selected[lead]={'status':'reference_match_count_'+str(len(hits))};continue
  item,bb,rows=hits[0];piece=int(item.filename.split(f'chr{chrom}.')[1].split('_')[0]);d=work/lead;(d/'bim').mkdir(parents=True);(d/'ld').mkdir();(d/'bim'/Path(item.filename).name).write_bytes(bb)
  selected[lead]={'status':'selected','chr':chrom,'piece':piece,'bim_member':item.filename,'bim_sha256':hashlib.sha256(bb).hexdigest(),'n_reference':len(rows),'min_pos':min(int(r[3]) for r in rows),'max_pos':max(int(r[3]) for r in rows),'snps':[r[1] for r in rows]}
size=1830679554;url=base+'LDfile.zip/content'
class Remote(io.RawIOBase):
 def __init__(self):self.pos=0
 def seekable(self):return True
 def readable(self):return True
 def tell(self):return self.pos
 def seek(self,n,whence=0):self.pos=n if whence==0 else self.pos+n if whence==1 else size+n;return self.pos
 def read(self,n=-1):
  if n<0:n=size-self.pos
  n=min(n,size-self.pos)
  if n<=0:return b''
  r=requests.get(url,headers={'Range':f'bytes={self.pos}-{self.pos+n-1}'},params={'start':self.pos},timeout=90);r.raise_for_status()
  assert r.status_code==206 and r.headers.get('Content-Range')==f'bytes {self.pos}-{self.pos+n-1}/{size}' and len(r.content)==n
  self.pos+=n;return r.content
with zipfile.ZipFile(Remote()) as z:
 for lead,v in selected.items():
  if v['status']!='selected':continue
  names=[n for n in z.namelist() if n.startswith(f"LD.path/ukb_chr{v['chr']}.{v['piece']}_") and n.endswith('_LDSVD.rda')]
  assert len(names)==1
  item=z.getinfo(names[0]);assert item.file_size < 8_000_000 and shutil.disk_usage(work).free>item.file_size+8_000_000
  data=z.read(item);(work/lead/'ld'/Path(names[0]).name).write_bytes(data)
  v.update(ld_member=names[0],ld_bytes=len(data),ld_crc32=f'{item.CRC:08x}',ld_sha256=hashlib.sha256(data).hexdigest())
  print(lead,v['chr'],v['piece'],len(data),flush=True)
inputbase=Path('/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/frailty_v1/lava_sensitivity_v1/inputs')
for lead,v in selected.items():
 if v['status']!='selected':continue
 for trait in ('frailty','insomnia'):
  snps=set(v['snps']);source=inputbase/trait/f"chr{v['chr']:02d}.sumstats.tsv.gz"
  with gzip.open(source,'rt') as src,(work/lead/f'{trait}.tsv').open('w') as dst:
   reader=csv.DictReader(src,delimiter='\t');writer=csv.DictWriter(dst,fieldnames=reader.fieldnames,delimiter='\t',lineterminator='\n');writer.writeheader();count=0
   for row in reader:
    if row['SNP'] in snps:writer.writerow(row);count+=1
  v[trait+'_rows']=count
 del v['snps']
(out/'prepare_receipt.json').write_text(json.dumps({'reference':'https://doi.org/10.5281/zenodo.14825987','bim_archive_md5':'608f037e4fc4d69f0ee9231c7fc005bb','LD_archive_fully_verified':False,'work':str(work),'blocks':selected},indent=2)+'\n')
