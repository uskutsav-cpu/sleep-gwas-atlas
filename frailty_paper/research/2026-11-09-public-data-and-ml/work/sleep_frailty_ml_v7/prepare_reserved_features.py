"""Physical final-study input only. No result table is opened."""
from common import *
from sequence_features import Encoder
def main():
    audit=json.loads((OUT/'reserved_design_and_original_source_audit.json').read_text());path=DATA/'reserved_physical_design.csv'
    rows=list(csv.DictReader(path.open()));x=Encoder().features([r['reference'] for r in rows],[r['alternate'] for r in rows],[int(r['changed_position']) for r in rows]);out=DATA/'reserved_feature_bank.npy';np.save(out,x)
    write_json(OUT/'reserved_feature_receipt.json',{'rows':len(rows),'columns':x.shape[1],'path':str(out),'sha256':sha(out),'design_sha256':sha(path),'feature_code_sha256':sha(WORK/'sequence_features.py'),'outcomes_read':False});print('RESERVED INPUT FEATURES READY',len(rows),flush=True)
if __name__=='__main__':main()
