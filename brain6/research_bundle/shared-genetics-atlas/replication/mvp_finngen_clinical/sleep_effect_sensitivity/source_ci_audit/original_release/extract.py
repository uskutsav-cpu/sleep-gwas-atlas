import csv,hashlib,io,json,re,zlib
from pathlib import Path
import requests
p=Path(__file__).resolve().parent
base='https://ftp.ncbi.nlm.nih.gov/dbgap/studies/phs002453/phs002453.v1.p1/analyses/GIA/phs002453.MVP_R4.1000G_AGR.GIA.PheCodes_Neurological.analysis-PI.MULTI.tar'
name='MVP_R4.1000G_AGR.GIA.PheCodes_Neurological/MVP_R4.1000G_AGR.Phe_327_4.EUR.GIA.dbGaP.txt.gz'
meta=name.replace('.txt.gz','.metadata.txt');resp=requests.get(base+'.table_of_contents.txt',timeout=45);resp.raise_for_status();toc=resp.text;offset=0;found={}
for line in toc.splitlines():
 m=re.match(r'\S+\s+\S+\s+(\d+)\s+\S+\s+\S+\s+(.+)$',line);assert m
 length=int(m.group(1));n=m.group(2)
 if n in (name,meta):found[n]=(offset,length)
 offset+=512+((length+511)//512)*512
assert set(found)=={name,meta} and found[name][1]==578841528
def getrange(start,end,stream=False):
 r=requests.get(base,headers={'Range':f'bytes={start}-{end}'},params={'start':start},timeout=(30,90),stream=stream);r.raise_for_status();assert r.status_code==206 and r.headers['Content-Range']==f'bytes {start}-{end}/180631818240';return r
def header(n):
 off,size=found[n];raw=getrange(off,off+511).content;assert len(raw)==512 and raw[:100].split(b'\0')[0].decode()==n
 ck=int(raw[148:156].strip(b' \0'),8);assert sum(raw[:148])+8*32+sum(raw[156:])==ck;assert int(raw[124:136].strip(b' \0'),8)==size
 return off+512,size
mstart,msize=header(meta);m=getrange(mstart,mstart+msize-1).content;assert len(m)==msize;(p/'source_metadata.txt').write_bytes(m)
start,size=header(name);r=getrange(start,start+size-1,True);hasher=hashlib.sha256();decomp=zlib.decompressobj(31);pending=b'';num=0;members=0;matches=[];targets={'rs2582894','rs7160630','rs62264767','rs2431108','rs6452785','rs1476535','rs62264764'};columns=None;downloaded=0;at_boundary=True
for chunk in r.iter_content(chunk_size=2**20):
 if not chunk:continue
 downloaded+=len(chunk);hasher.update(chunk);remaining=chunk
 while remaining:
  decoded=decomp.decompress(remaining);remaining=decomp.unused_data;at_boundary=decomp.eof
  data=pending+decoded;lines=data.split(b'\n');pending=lines.pop()
  for line in lines:
   if columns is None:columns=line.decode().split('\t');continue
   num+=1
   if line.split(b'\t',1)[0].decode() in targets:
    row=dict(zip(columns,line.decode().split('\t')));matches.append(row);print('FOUND',row['SNP_ID'],flush=True)
  if decomp.eof:
   members+=1;decomp=zlib.decompressobj(31)
  else:break
 if downloaded % (64*2**20)==0:print('STREAMED',downloaded,'bytes',flush=True)
assert downloaded==size and members>0 and at_boundary and not pending
(p/'original_rows.json').write_text(json.dumps(matches,indent=2)+'\n')
(p/'receipt.json').write_text(json.dumps({'archive':base,'tar_member':name,'tar_header_offset':start-512,'member_bytes':size,'member_sha256':hasher.hexdigest(),'table_of_contents_sha256':hashlib.sha256(resp.content).hexdigest(),'metadata_sha256':hashlib.sha256(m).hexdigest(),'rows':num,'targets':sorted(targets),'matches':[x['SNP_ID'] for x in matches],'bgzf_members_crc_checked':members,'full_tar_verified':False},indent=2)+'\n');print('COMPLETE',downloaded,num,flush=True)
