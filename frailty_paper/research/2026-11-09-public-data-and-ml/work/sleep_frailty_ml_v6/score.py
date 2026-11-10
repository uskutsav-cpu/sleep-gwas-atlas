"""Score paired DNA sequences, preserving explicit uncertainty/abstention.

Example:
PYTHONPATH=work/stage4_env/lib/python3.13/site-packages python3 \
 work/sleep_frailty_ml_v6/score.py --input pairs.csv --output scores.csv
Input columns: id, reference, alternate, cell_group. Optional genome_build
must be GRCh38/hg38. Exactly one substitution, 151–1024 A/C/G/T bases.
This predicts a laboratory statistical detection; never patient disease risk.
"""
import argparse
from common import *
from encode_inputs import features
from predict_checkpoints import predict_registry
from calibrate_uncertainty import apply_policy
def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--input',required=True);parser.add_argument('--output',required=True);args=parser.parse_args();registry_path=OUT/'final_fit_registry.json';policy=json.loads((OUT/'uncertainty_policy.json').read_text());assert sha(registry_path)==policy['primary_registry_sha256'];registry=json.loads(registry_path.read_text());rows=list(csv.DictReader(Path(args.input).open()));valid=[];output=[]
    for i,r in enumerate(rows):
        a=r.get('reference','').upper();b=r.get('alternate','').upper();changes=[j for j,(x,y) in enumerate(zip(a,b)) if x!=y]
        base={'id':r.get('id',str(i)),'cell_group':r.get('cell_group',''),'status':'unranked_invalid_input','ranked_score':'','calibrated_laboratory_call_probability':'','prediction_set':''}
        if len(a)==len(b) and 151<=len(a)<=1024 and set(a+b)<=set('ACGT') and len(changes)==1 and r.get('genome_build','GRCh38') in ['GRCh38','hg38']:
            valid.append((i,a,b,changes[0],base['cell_group']))
        output.append(base)
    if valid:
        x=features([r[1] for r in valid],[r[2] for r in valid],[r[3] for r in valid]);cells=np.array([r[4] for r in valid]);prediction,members,distance=predict_registry(registry,x,np.arange(len(valid)),cells);result=apply_policy(prediction['primary'],cells,distance,policy)
        for j,r in enumerate(valid):
            i=r[0];output[i]['status']=str(result['status'][j]);output[i]['ranked_score']=float(result['ranked_score'][j]) if np.isfinite(result['ranked_score'][j]) else '';output[i]['calibrated_laboratory_call_probability']=float(result['calibrated_call_probability'][j]);output[i]['prediction_set']=','.join(str(k) for k in [0,1] if result['set_includes'+str(k)][j])
    with Path(args.output).open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(output[0]) if output else ['id','cell_group','status','ranked_score','calibrated_laboratory_call_probability','prediction_set']);w.writeheader();w.writerows(output)
    print('Scored',len(rows),'pairs;',sum(r['status']=='ranked' for r in output),'ranked; others explicitly unranked.')
if __name__=='__main__':main()
