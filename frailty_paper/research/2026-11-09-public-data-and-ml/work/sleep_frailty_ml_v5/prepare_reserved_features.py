"""Encode the reserved study's audited design without opening measurements."""
from common import *
from prepare_blind_test_design import fasta
from encode_inputs import features
def main():
    rows=list(csv.DictReader((DATA/'kosicki_blind_design.csv').open()));design_audit=json.loads((OUT/'kosicki_blind_design_audit.json').read_text());assert sha(DATA/'kosicki_blind_design.csv')==design_audit['design_sha256'];target=DATA/'kosicki_fixed_features.npy';receipt_path=OUT/'kosicki_fixed_features_receipt.json'
    if target.exists() and receipt_path.exists():
        receipt=json.loads(receipt_path.read_text());assert sha(target)==receipt['sha256'];print('Reused byte-verified reserved design features');return
    seq=fasta(SOURCES/'kosicki_ENCFF752NFF.fasta.gz');ref=[seq[r['reference_name']] for r in rows];alt=[seq[r['name']] for r in rows];pos=[int(r['changed_position']) for r in rows];x=features(ref,alt,pos);np.save(target,x)
    write_json(receipt_path,{'file':str(target),'sha256':sha(target),'shape':list(x.shape),'design_sha256':sha(DATA/'kosicki_blind_design.csv'),'source_feature_receipt_sha256':sha(OUT/'fixed_features_receipt.json'),'source_base_sha256':sha(OLD/'hyenadna/model.safetensors'),'measured_labels_opened_or_used':False,'scope':'Fixed encoding of the entire designed final-study sequence universe, not a biological evaluation or model fit.'});print('RESERVED DESIGN FEATURES COMPLETE',x.shape,flush=True)
if __name__=='__main__':main()
