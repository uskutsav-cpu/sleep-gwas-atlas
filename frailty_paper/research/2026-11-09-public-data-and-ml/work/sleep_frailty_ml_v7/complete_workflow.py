"""Sequential resumable execution; stops on a failed scientific/implementation gate."""
import argparse,subprocess,sys
from common import *
def main():
    p=argparse.ArgumentParser();p.add_argument('--wait-for-training-pid',type=int);args=p.parse_args()
    if args.wait_for_training_pid:
        while subprocess.run(['ps','-p',str(args.wait_for_training_pid)],capture_output=True).returncode==0:time.sleep(10)
    expected=['raw_MAE','raw_MAE_pairs','raw_MAE_no_cell_pairs']
    for recipe in expected:
        for seed in [42,73,101]:assert (OUT/f'development_{recipe}_{seed}_record.json').exists(),f'Training incomplete:{recipe}/{seed}'
    for name in ['finalize.py','audit_partitions.py','evaluate.py','audit_models.py','integrate.py','build_delivery.py','audit_completion.py']:
        if name=='audit_partitions.py' and (OUT/'independent_partition_audit.json').exists():
            audit=json.loads((OUT/'independent_partition_audit.json').read_text())
            if audit['passed'] and all(sha(DATA/n)==h for n,h in audit.get('input_sha256',{}).items()) and len(audit.get('input_sha256',{}))==3:continue
        print('RUNNING',name,time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),flush=True);subprocess.run([sys.executable,str(WORK/name)],env=environment(),check=True)
    print('MODELLING WORKFLOW FINISHED; REVIEW SCIENTIFIC REQUIREMENT AUDIT',flush=True)
if __name__=='__main__':main()
