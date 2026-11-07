#!/usr/bin/env python3
"""Submit one permitted top-AVI variant per remaining region to JBrowse.

Keep status JSON under scratch because it includes short-lived signed URLs.
If submission outcome is unknown, stop for operator review rather than re-submit.
"""
import argparse,csv,json,ssl,subprocess,sys,time,urllib.request
from pathlib import Path
P=Path(__file__).resolve().parent; BASE='https://0t0e9nn6bj.execute-api.us-east-2.amazonaws.com/Prod/api/predict'; CTX=ssl.create_default_context(cafile='/etc/ssl/cert.pem')
p=argparse.ArgumentParser();p.add_argument('--scratch',type=Path,required=True);p.add_argument('--output-root',type=Path,required=True);args=p.parse_args();args.scratch.mkdir(parents=True,exist_ok=True)
with (P/'results/avi_mirror/variant_evidence_avi.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
groups={}
for r in rows:groups.setdefault(r['candidate_locus_id'],[]).append(r)
if len(groups)!=6:raise ValueError('Expected six regions')
for cid,group in sorted(groups.items()):
 if cid=='placo_insomnia__adhd_chr5_103447968_104447968':continue
 r=sorted(group,key=lambda x:(-int(x['avi_raw_e5']),x['atlas_variant']))[0];dirname=args.scratch/cid;dirname.mkdir(exist_ok=True); submission=dirname/'submission.json';status=dirname/'status.json'; output=args.output_root/cid
 if output.exists():print(cid,'already extracted',flush=True);continue
 pos=int(r['pos_grch38_1based']);req={'mode':'variant','region':{'chromosome':r['chromosome'],'start':pos-1-5000,'end':pos-1+5000,'strand':'.'},'sequence_length':1048576,'output_types':['rna_seq','dnase'],'ontology_terms':['UBERON:0001870'],'organism':'human','variant':{'chromosome':r['chromosome'],'position':pos,'reference_bases':r['ref_grch38'],'alternate_bases':r['alt_grch38'],'name':r['snp']}}
 if not submission.exists():
  (dirname/'request.json').write_text(json.dumps(req,indent=2)+'\n')
  # Marker prevents automatic duplicate submissions after an unknown outcome.
  marker=dirname/'SUBMISSION_IN_PROGRESS'
  if marker.exists():raise RuntimeError(f'Unknown prior submission outcome for {cid}; investigate before retrying')
  marker.write_text('A request was initiated. Do not retry automatically.\n')
  data=json.dumps(req).encode();http=urllib.request.Request(BASE,data=data,headers={'Content-Type':'application/json'},method='POST')
  with urllib.request.urlopen(http,context=CTX,timeout=40) as response:result=json.load(response)
  if not result.get('token'):raise ValueError('Submission returned no token')
  submission.write_text(json.dumps(result,indent=2)+'\n');marker.unlink()
  print('submitted',cid,r['snp'],flush=True)
 result=json.loads(submission.read_text());token=result['token']
 for attempt in range(90):
  with urllib.request.urlopen(BASE+'/'+token,context=CTX,timeout=25) as response:state=json.load(response)
  status.write_text(json.dumps(state)+'\n')
  if state.get('status')=='ready' and state.get('manifest'):break
  if state.get('status')=='error':raise RuntimeError(f'Service reported error for {cid}')
  time.sleep(2)
 else:raise TimeoutError(f'Prediction still pending for {cid}; rerun resumes same token')
 print('ready',cid,flush=True)
 subprocess.run([sys.executable,str(P/'extract_jbrowse_window.py'),'--status',str(status),'--out',str(output),'--expected-variant',r['atlas_variant']],check=True,stdout=subprocess.DEVNULL)
 print('extracted',cid,flush=True)
