#!/usr/bin/env python3
"""Bounded public EBI v5b metadata and exact HTTP-range acquisition."""
import argparse, datetime, gzip, hashlib, json, pathlib, re, shutil, struct, subprocess
HERE=pathlib.Path(__file__).resolve().parent
ROOT=HERE.parents[3]
RAW=ROOT/'work'/'ld_genotypes_research_v1'
BASE='https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/'
OLD=pathlib.Path('/Users/swethasunilkumar/Documents/Codex/2026-09-19/can/work/sleep-gwas-atlas-v03')
REGIONS=[(5,103447968,104447968),(6,100630147,102636772),(11,112459489,114257728)]
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def dump(path,obj):path.write_text(json.dumps(obj,indent=2,sort_keys=True)+'\n')
def fetch(url,path,cap,byte_range=None):
 path.parent.mkdir(parents=True,exist_ok=True)
 headers=path.with_suffix(path.suffix+'.headers')
 cmd=['curl','--fail','--silent','--show-error','--location','--max-time','60','--max-filesize',str(cap),'-D',str(headers),'-o',str(path)]
 if byte_range:cmd+=['--range',f'{byte_range[0]}-{byte_range[1]}']
 result=subprocess.run(cmd+[url],capture_output=True,text=True)
 rec={'url':url,'path':str(path),'retrieved_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'returncode':result.returncode,'stderr':result.stderr,'range':byte_range,'headers_sha256':sha(headers) if headers.exists() else None}
 if path.exists():rec.update(bytes=path.stat().st_size,sha256=sha(path))
 rec['response_headers']=headers.read_text() if headers.exists() else ''
 if result.returncode!=0:raise RuntimeError(json.dumps(rec))
 if path.stat().st_size>cap:raise RuntimeError('Exceeded cap')
 if byte_range:
  expected=f'bytes {byte_range[0]}-{byte_range[1]}/'
  if not re.search(r'HTTP/\S+ 206\b',rec['response_headers']) or expected not in rec['response_headers'] or path.stat().st_size!=byte_range[1]-byte_range[0]+1:raise RuntimeError('Rejected nonmatching range response '+json.dumps(rec))
 return rec

def index_ranges(path,chrom,start,stop):
 d=gzip.decompress(path.read_bytes());o=4
 if d[:4]!=b'TBI\1':raise ValueError('Invalid TBI')
 h=struct.unpack_from('<8i',d,o);o+=32
 if h[0]!=1 or h[1]&0xffff!=2:raise ValueError('Expected one-reference VCF TBI')
 names=d[o:o+h[7]].rstrip(b'\0').decode().split('\0');o+=h[7]
 bins={0};end=stop-1
 for base,shift in [(1,26),(9,23),(73,20),(585,17),(4681,14)]:bins.update(range(base+((start-1)>>shift),base+(end>>shift)+1))
 nb=struct.unpack_from('<i',d,o)[0];o+=4;chunks=[]
 for _ in range(nb):
  binid,n=struct.unpack_from('<Ii',d,o);o+=8
  values=[struct.unpack_from('<QQ',d,o+16*j) for j in range(n)];o+=16*n
  if binid in bins:chunks.extend(values)
 nl=struct.unpack_from('<i',d,o)[0];o+=4;lin=struct.unpack_from(f'<{nl}Q',d,o)
 if names!=[str(chrom)]:raise ValueError('TBI chromosome mismatch')
 minimum=lin[(start-1)>>14]
 pairs=sorted((max(s,minimum)>>16,(t>>16)+65535) for s,t in chunks if t>minimum)
 merged=[]
 for a,b in pairs:
  if merged and a<=merged[-1][1]+1:merged[-1][1]=max(b,merged[-1][1])
  else:merged.append([a,b])
 return merged

def metadata():
 RAW.mkdir(parents=True,exist_ok=True)
 files=['README_phase3_callset_20150220','README_known_issues_20200731','README_integrated_call_samples_v3.20250704.ALL','integrated_call_samples_v3.20130502.ALL.panel','integrated_call_samples_v3.20250704.ALL.ped']
 files += [f'ALL.chr{c}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi' for c,a,b in REGIONS]
 manifest={'protocol_sha256':sha(HERE/'LD_PRE_OUTCOME_PROTOCOL.md'),'frozen_before_bodies':True,'total_metadata_cap':1000000,'sources':[{'url':BASE+x,'output':str(RAW/x),'cap':300000} for x in files]}
 dump(HERE/'metadata_acquisition_manifest.json',manifest)
 receipts=[]
 for obj in manifest['sources']:
  receipts.append(fetch(obj['url'],pathlib.Path(obj['output']),obj['cap']))
  dump(HERE/'metadata_acquisition_receipts.json',receipts)
  if sum(r['bytes'] for r in receipts)>manifest['total_metadata_cap']:raise RuntimeError('Metadata total cap exceeded')
 for c,a,b in REGIONS:
  name=f'ALL.chr{c}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi'
  if sha(RAW/name)!=sha(OLD/name):raise RuntimeError('Retained/current TBI mismatch: '+name)
 plan={'protocol_sha256':manifest['protocol_sha256'],'metadata_receipts_sha256':sha(HERE/'metadata_acquisition_receipts.json'),'frozen_before_genotype_bodies':True,'max_compressed_bytes':30000000,'regions':[]}
 for c,a,b in REGIONS:
  name=f'ALL.chr{c}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz'
  ranges=index_ranges(RAW/(name+'.tbi'),c,a,b)
  plan['regions'].append({'chr':c,'start':a,'stop':b,'url':BASE+name,'index_sha256':sha(RAW/(name+'.tbi')),'header_range':[0,65535],'genotype_ranges':ranges})
 plan['total_compressed_bytes']=sum(65536+sum(y-x+1 for x,y in r['genotype_ranges']) for r in plan['regions'])
 if plan['total_compressed_bytes']>plan['max_compressed_bytes']:raise RuntimeError('Genotype acquisition total cap exceeded')
 dump(HERE/'genotype_acquisition_manifest.json',plan)
 print(json.dumps({'metadata_bytes':sum(r['bytes'] for r in receipts),'planned_genotype_bytes':plan['total_compressed_bytes'],'protocol_sha256':plan['protocol_sha256']}))

def download():
 plan=json.loads((HERE/'genotype_acquisition_manifest.json').read_text())
 if sha(HERE/'LD_PRE_OUTCOME_PROTOCOL.md')!=plan['protocol_sha256']:raise RuntimeError('Protocol changed after acquisition freeze')
 if shutil.disk_usage(RAW).free<2**30:raise RuntimeError('Less than1GiB free')
 receipts=[]
 for obj in plan['regions']:
  for typ,ranges in [('header',[obj['header_range']]),('genotypes',obj['genotype_ranges'])]:
   for i,r in enumerate(ranges):
    path=RAW/f'chr{obj["chr"]}_{typ}_{i}_{r[0]}_{r[1]}.bgzf'
    receipts.append(fetch(obj['url'],path,r[1]-r[0]+1,r))
    dump(HERE/'genotype_acquisition_receipts.json',receipts)
 print(json.dumps({'retrieved_ranges':len(receipts),'bytes':sum(r['bytes'] for r in receipts)}))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['metadata','download']);a=p.parse_args();globals()[a.stage]()
