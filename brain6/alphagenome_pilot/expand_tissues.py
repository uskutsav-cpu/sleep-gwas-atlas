#!/usr/bin/env python3
"""Resumable sequential 1 Mb AlphaGenome frontal-cortex pilot extension.

Each 2,001-bp pair of float32 tracks is stored as gzip, not full model arrays.
Scratch contains expiring signed URLs and is never copied into deliverables.
A pending marker protects against repeating a POST with an unknown outcome.
"""
import argparse,csv,gzip,hashlib,json,math,os,shutil,ssl,struct,sys,time,urllib.error,urllib.parse,urllib.request
from pathlib import Path
from datetime import datetime,timezone
P=Path(__file__).resolve().parent; E=P/'results/avi_mirror/variant_evidence_avi.tsv';BASE='https://0t0e9nn6bj.execute-api.us-east-2.amazonaws.com/Prod/api/predict';CTX=ssl.create_default_context(cafile='/etc/ssl/cert.pem'); N=2001

def sha(b):return hashlib.sha256(b).hexdigest()
def jread(p):return json.loads(p.read_text())
def atomic(p,data):
 tmp=p.with_name(p.name+'.tmp');tmp.write_bytes(data);os.replace(tmp,p)
def jwrite(p,obj):atomic(p,(json.dumps(obj,sort_keys=True,indent=2)+'\n').encode())
def http_json(url,data=None):
 req=urllib.request.Request(url,data=data,headers={'Content-Type':'application/json'} if data else {})
 with urllib.request.urlopen(req,context=CTX,timeout=45) as r:return json.load(r)
def variant_string(v):return f"{v['chromosome']}:{v['position']}:{v['reference_bases']}>{v['alternate_bases']}"
def extract(state,expected):
 m=state['manifest'];v=m['variant']; iv=m['interval'];p0=v['position']-1
 if state['status']!='ready' or m['mode']!='variant' or m['organism']!='human' or m['sequence_length']!=1048576 or m['ontology_terms']!=['UBERON:0001870'] or variant_string(v)!=expected:raise ValueError('prediction identity, model context or tissue mismatch')
 offset=p0-1000-iv['start'];
 if offset<0 or p0+1001>iv['end']:raise ValueError('window outside prediction')
 result=[]; data=[]; hashes=[]
 for name in ('dnase','rna_seq'):
  ref=m['outputs'][name]['reference'];alt=m['outputs'][name]['alternate']
  if ref['resolution']!=1 or alt['resolution']!=1 or ref['num_bins']!=1048576 or ref['tracks']!=alt['tracks']:raise ValueError('track metadata mismatch')
  for i,track in enumerate(ref['tracks']):
   pairs=[]
   for mode,entry in [('reference',ref),('alternate',alt)]:
    start=i*entry['bytes_per_track']+offset*4;last=start+N*4-1
    req=urllib.request.Request(entry['url'],headers={'Range':f'bytes={start}-{last}'})
    with urllib.request.urlopen(req,context=CTX,timeout=45) as response:
     body=response.read();status=response.status;range_header=response.headers.get('Content-Range','')
    if status!=206 or not range_header.startswith(f'bytes {start}-{last}/') or len(body)!=N*4:raise ValueError('invalid byte range')
    values=struct.unpack('<'+str(N)+'f',body)
    if not all(math.isfinite(x) for x in values):raise ValueError('non-finite prediction')
    pairs.append(values);hashes.append({'output_type':name,'track_ordinal':i,'variant_type':mode,'byte_range':f'{start}-{last}','sha256':sha(body)});data.append(body)
   r,b=pairs;sr=sum(r[750:1251]);sb=sum(b[750:1251]);ratio=math.log2((sb+1)/(sr+1)) if sr+1>0 and sb+1>0 else None
   meta={k:x for k,x in track.items() if k not in ('url','key')}
   result.append({'output_type':name,'track_ordinal':i,'metadata':meta,'ref_sum_501bp':sr,'alt_sum_501bp':sb,'local_log2_ratio_plus_one':ratio})
 if len(result)!=4:raise ValueError('expected four output tracks')
 return {'variant':v,'interval':iv,'request':state['request'],'model_release':'not exposed','created':m['created'],'window_bp':[-1000,1000],'tracks':result,'source_ranges':hashes},b''.join(data)
def save(output,meta,data):
 output.mkdir(parents=True,exist_ok=True)
 packed=gzip.compress(data,compresslevel=6,mtime=0);atomic(output/'track_window.f32.gz',packed);meta['window_format']='little-endian float32, tracks in metadata order, reference then alternate, 2001 values per array';meta['track_window_sha256']=sha(packed);jwrite(output/'provenance.json',meta)
def import_existing(row,folder,output):
 meta=jread(folder/'provenance.json');v=meta['request']['variant']
 if variant_string(v)!=row['atlas_variant']:raise ValueError('existing variant mismatch')
 for name,h in meta['output_sha256'].items():
  if sha((folder/name).read_bytes())!=h:raise ValueError('existing prediction hash mismatch')
 with (folder/'track_window.tsv').open() as f:window=list(csv.DictReader(f,delimiter='\t'))
 with (folder/'track_summary.tsv').open() as f:summary=list(csv.DictReader(f,delimiter='\t'))
 chunks=[]
 for t in meta['tracks']:
  part=[x for x in window if x['output_type']==t['output_type'] and x['track_ordinal']==str(t['track_ordinal'])]
  if len(part)!=N or [int(x['offset_from_variant_bp']) for x in part]!=list(range(-1000,1001)):raise ValueError('existing window mismatch')
  for mode in ('reference','alternate'):chunks.append(struct.pack('<'+str(N)+'f',*(float(x[mode]) for x in part)))
 meta2={'variant':v,'request':meta['request'],'interval':meta['request']['interval'],'model_release':'not exposed','created':meta['manifest_created_epoch'],'window_bp':[-1000,1000],'tracks':[], 'source_ranges':meta['range_sha256'],'imported_from':str(folder.relative_to(P)),'imported_provenance_sha256':sha((folder/'provenance.json').read_bytes())}
 meta2['candidate_locus_id']=row['candidate_locus_id'];meta2['evidence_sha256']=sha(E.read_bytes())
 for t,s in zip(meta['tracks'],summary):
  meta2['tracks'].append({'output_type':t['output_type'],'track_ordinal':t['track_ordinal'],'metadata':t,'ref_sum_501bp':float(s['ref_sum_501bp']),'alt_sum_501bp':float(s['alt_sum_501bp']),'local_log2_ratio_plus_one':float(s['local_log2_ratio_plus_one']) if s['local_log2_ratio_plus_one'] else None})
 save(output,meta2,b''.join(chunks))
def main():
 p=argparse.ArgumentParser();p.add_argument('--scratch',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--max-new',type=int,default=425);p.add_argument('--delay',type=float,default=1.0);p.add_argument('--region');a=p.parse_args();a.scratch.mkdir(parents=True,exist_ok=True);a.output.mkdir(parents=True,exist_ok=True)
 rows=list(csv.DictReader(E.open(),delimiter='\t'))
 if len(rows)!=425 or len({r['atlas_variant'] for r in rows})!=425 or any(r['mapping_status']!='READY_UNIQUE_REFERENCE_VERIFIED' for r in rows):raise ValueError('input universe mismatch')
 rows.sort(key=lambda r:(not r['candidate_locus_id'].startswith('placo_insomnia__adhd_chr5_103447968'),r['candidate_locus_id'],r['pos_grch38_1based'],r['snp']))
 fhash=sha(E.read_bytes());binding=a.output/'input_binding.json'
 if binding.exists() and jread(binding)['evidence_sha256']!=fhash:raise ValueError('input changed')
 if not binding.exists():jwrite(binding,{'evidence_sha256':fhash,'expected_variant_count':425,'ontology':'UBERON:0001870','outputs':['dnase','rna_seq'],'sequence_length':1048576,'selection':'all 425 prepared credible-set union variants','noncommercial':True})
 if a.region:
  rows=[r for r in rows if r['candidate_locus_id']==a.region]
  if not rows:raise ValueError('Requested region not in input universe')
 existing={}
 for folder in [P/'results/jbrowse_chr5',*(P/'results/jbrowse_five').iterdir()]:
  meta=jread(folder/'provenance.json');existing[variant_string(meta['request']['variant'])]=folder
 new=0
 for i,row in enumerate(rows):
  name=sha(row['atlas_variant'].encode())[:20];out=a.output/'variants'/name;state_dir=a.scratch/name;state_dir.mkdir(parents=True,exist_ok=True)
  if (out/'provenance.json').exists():
   m=jread(out/'provenance.json')
   if variant_string(m['variant'])!=row['atlas_variant'] or sha((out/'track_window.f32.gz').read_bytes())!=m['track_window_sha256']:raise ValueError('existing output invalid')
   continue
  if shutil.disk_usage(a.output).free<80*1024*1024:print('STOP_LOW_DISK',flush=True);return 2
  if row['atlas_variant'] in existing:import_existing(row,existing[row['atlas_variant']],out);print(f'{i+1}/425 IMPORTED {row["snp"]}',flush=True);continue
  if new>=a.max_new:print('STOP_BATCH_LIMIT',flush=True);return 0
  submission=state_dir/'submission.json';marker=state_dir/'SUBMISSION_UNKNOWN'
  pos=int(row['pos_grch38_1based']);req={'mode':'variant','region':{'chromosome':row['chromosome'],'start':pos-5001,'end':pos+4999,'strand':'.'},'sequence_length':1048576,'output_types':['rna_seq','dnase'],'ontology_terms':['UBERON:0001870'],'organism':'human','variant':{'chromosome':row['chromosome'],'position':pos,'reference_bases':row['ref_grch38'],'alternate_bases':row['alt_grch38'],'name':row['snp']}}
  if not submission.exists():
   if marker.exists():print(f'STOP_UNKNOWN_SUBMISSION {row["snp"]}',flush=True);return 3
   jwrite(state_dir/'request.json',req);marker.write_text('Submission started; outcome must be resolved before retry.\n')
   try: answer=http_json(BASE,json.dumps(req).encode())
   except urllib.error.HTTPError as exc:
    # 4xx indicates rejection; 5xx can have an unknown submission outcome.
    jwrite(state_dir/'http_error.json',{'code':exc.code,'when':datetime.now(timezone.utc).isoformat()})
    if 400<=exc.code<500:marker.unlink()
    print(f'STOP_HTTP_{exc.code} {row["snp"]}',flush=True);return 4
   except Exception as exc:print(f'STOP_UNKNOWN_SUBMISSION {row["snp"]} {type(exc).__name__}',flush=True);return 3
   if not answer.get('token'):print('STOP_NO_TOKEN',flush=True);return 3
   jwrite(submission,answer);marker.unlink();new+=1
  token=jread(submission)['token']
  for attempt in range(90):
   try: state=http_json(BASE+'/'+urllib.parse.quote(token,safe=''))
   except Exception as exc:print(f'STOP_POLL_RETRY_SAME_TOKEN {row["snp"]} {type(exc).__name__}',flush=True);return 5
   if state.get('status')=='ready' and state.get('manifest'):break
   if state.get('status')=='error':jwrite(state_dir/'terminal_status.json',{'status':'error'});print(f'STOP_SERVICE_ERROR {row["snp"]}',flush=True);return 6
   time.sleep(2)
  else:print(f'STOP_POLL_TIMEOUT_SAME_TOKEN {row["snp"]}',flush=True);return 5
  try:meta,data=extract(state,row['atlas_variant']);meta['candidate_locus_id']=row['candidate_locus_id'];meta['evidence_sha256']=fhash;save(out,meta,data)
  except Exception as exc:print(f'STOP_EXTRACT_RETRY_SAME_TOKEN {row["snp"]} {type(exc).__name__}',flush=True);return 7
  print(f'{i+1}/425 SCORED {row["snp"]}',flush=True);time.sleep(a.delay)
 return 0
if __name__=='__main__':sys.exit(main())
