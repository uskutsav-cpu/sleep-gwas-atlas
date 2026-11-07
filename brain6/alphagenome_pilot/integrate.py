#!/usr/bin/env python3
"""Join a completed Atlas run to Brain6; preserve failed/unattempted statuses."""
import argparse,csv,gzip,hashlib,json
from pathlib import Path

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def read(path):
    with path.open(newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--evidence',type=Path,required=True);p.add_argument('--run',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    cfg=json.loads((a.run/'run_config.json').read_text())
    if sha(a.evidence)!=cfg['evidence_sha256']:raise ValueError('Evidence file differs from scored input')
    statuses=read(a.run/'query_status.tsv'); states={s['variant']:s for s in statuses}
    if len(states)!=len(statuses):raise ValueError('Duplicate variant status')
    scored={}
    for v,s in states.items():
        if s['status'] not in {'SCORED','NO_RETURNED_TRACKS','QUERY_FAILED'}:raise ValueError('Unknown score status')
        if s['status'] in {'SCORED','NO_RETURNED_TRACKS'}:
            path=(a.run/s['file']).resolve()
            if path.parent!=(a.run/'scores').resolve():raise ValueError('Invalid result path')
            with gzip.open(path,'rt',encoding='utf-8') as f: payload=json.load(f)
            if payload['variant']!=v or len(payload['scores'])!=int(s['n_score_rows']):raise ValueError('Score result mismatch')
            if (s['status']=='SCORED') != bool(payload['scores']):raise ValueError('Score status inconsistent')
            scored[v]=payload['scores']
    evidence=read(a.evidence);base=list(evidence[0]);extras=['atlas_query_status','scorer','observation_metadata_json','track_metadata_json','raw_score','quantile_score']
    a.out.mkdir(parents=True,exist_ok=False)
    count=0
    with (a.out/'variant_molecular_evidence.tsv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=base+extras,delimiter='\t',lineterminator='\n');w.writeheader()
        for row in evidence:
            v=row['atlas_variant']; state=states.get(v)
            status=state['status'] if state else ('NOT_QUERIED' if v else 'MAPPING_HOLD')
            items=scored.get(v,[]) or [None]
            for score in items:
                extra={'atlas_query_status':status,'scorer':score['scorer'] if score else '', 'observation_metadata_json':json.dumps(score['observation'],sort_keys=True) if score else '', 'track_metadata_json':json.dumps(score['track'],sort_keys=True) if score else '', 'raw_score':score['raw_score'] if score and score['raw_score'] is not None else '', 'quantile_score':score['quantile_score'] if score and score['quantile_score'] is not None else ''}
                # The original field describes preparation time and must not be mistaken for this run's status.
                w.writerow({**row,**extra});count+=1
    receipt={'evidence_sha256':sha(a.evidence),'run_config_sha256':sha(a.run/'run_config.json'),'query_status_sha256':sha(a.run/'query_status.tsv'),'score_files_sha256':{s['file']:sha(a.run/s['file']) for s in statuses if s['file']},'output_sha256':sha(a.out/'variant_molecular_evidence.tsv'),'rows':count,'interpretation':'Predictions are exploratory molecular effects, not a new causal or replication classification.'}
    (a.out/'integration_receipt.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
    print(f'Wrote {count} evidence rows, including unqueried and held variants.')
if __name__=='__main__':main()
