#!/usr/bin/env python3
"""Explicit opt-in replay of bounded, public read-only acquisition manifests."""
import argparse,csv,hashlib,json,subprocess,time
from pathlib import Path
from audit_evidence import ROOT,AREA
p=argparse.ArgumentParser();p.add_argument('--execute',action='store_true');p.add_argument('--manifest',type=Path,required=True);a=p.parse_args()
items=list(csv.DictReader(a.manifest.open(),delimiter='\t'))
if not a.execute:print(json.dumps({'status':'DRY_RUN','objects':len(items),'manifest':str(a.manifest)}));raise SystemExit(0)
store=ROOT/'work/reacquisition';store.mkdir(parents=True,exist_ok=True);receipts=[]
for x in items:
 dest=store/(x['id']+'.object');cap=int(x.get('max_bytes') or 30000000)
 r=subprocess.run(['curl','--fail','--location','--silent','--show-error','--max-time','30','--max-filesize',str(cap),'--output',str(dest),x['url']],capture_output=True,text=True)
 d={**x,'timestamp_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'status':'RETRIEVED' if r.returncode==0 else 'UNAVAILABLE','stderr':r.stderr,'path':str(dest)}
 if not r.returncode:d.update(sha256=hashlib.sha256(dest.read_bytes()).hexdigest(),bytes=dest.stat().st_size)
 receipts.append(d);(store/'receipts.json').write_text(json.dumps(receipts,indent=2)+'\n')
