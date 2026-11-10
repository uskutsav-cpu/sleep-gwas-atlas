import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'sleep_frailty_ml_representation'))
from representation_common import *
REPRESENTATION=DEST
DEST=ROOT/'outputs/sleep_frailty_ml_native'
PREP=ROOT/'outputs/sleep_frailty_ai_v9_validation/next_source_screen/candidate_pretrained_readers'


def assert_versions_native():
    assert_versions_all()
    for name,h in json.loads((REPRESENTATION/'final_progress_verification.json').read_text())['files'].items():
        assert sha(REPRESENTATION/name)==h,name


def add_all_documented_signed_targets(d):
    first=load_signed_sidecar(d)
    folder=REPRESENTATION/'neural_source_trace';receipt=json.loads((folder/'neural_original_direction_recovery_receipt.json').read_text())
    assert receipt['passed'] and sha(folder/'neural_original_signed_development.npy')==receipt['sidecar_sha256']
    check=list(csv.DictReader((folder/'neural_direction_independent_R_receipt.tsv').open(),delimiter='\t'))[0]
    assert check['passed']=='TRUE' and int(check['all_recovered_rows_checked'])==receipt['recovered_original_signed_targets']
    values=np.load(folder/'neural_original_signed_development.npy');known=np.isfinite(values)
    assert not (known&np.isfinite(d['signed_effect'])).any() and not known[d['split']==2].any()
    assert np.allclose(abs(values[known]),d['effect'][known],rtol=1e-5,atol=1e-6)
    d['signed_effect'][known]=values[known]
    return first,known


def get_plan():
    plan={'recipes':['frozen_reader','adapt_last_stage'],'seed':42,'epochs':3,
      'phase_one_roles':['none','Salomon','Kosicki2025','Retallick2026'],
      'phase_two_roles_if_gate_passes':['Siraj','IGVF240K','Koesterich2023','McAfee2023'],
      'data_sha256':sha(OLD/'coverage_training_arrays.npz'),
      'documented_Siraj_signed_sidecar_sha256':sha(DIRECTION/'original_paired_signed_development.npy'),
      'documented_neuronal_signed_sidecar_sha256':sha(REPRESENTATION/'neural_source_trace/neural_original_signed_development.npy'),
      'pretrained_reader':'Original WTC11 author test1/validation2 checkpoint with audited published pretraining exposure; fixed K562/HepG2 features retained.',
      'reader_sha256':sha(PREP/'WTC11_author_test1_val2.safetensors'),
      'sequence_policy':'Complete actual source DNA at its real length and assay orientation. Native AGCT tokens, no invented flanks or padding. Only one documented changed base reconstructs the exact source ALT.',
      'magnitude_features':'42 symmetric scalar summaries plus actual length and changed position; 192 fixed projected global/local means and absolute allele differences. Native WTC11 summaries replace its original bank entries. Magnitude is invariant to exchanging alleles.',
      'model':'Same small CellMatchedEffectHead in both arms, 236 features; the original reader mapper/head and BN running statistics stay fixed.',
      'only_between_arm_difference':'Adaptation permits gradients/updates in reader main.3, with relative squared-distance anchor 0.01 to original stage weights. Frozen reader permits none. Every other reader parameter stays fixed.',
      'supervision':'Both arms use the same original magnitude, documented Kosicki/Siraj/Salomon signed coefficients, original allele activities and nearby/cross-cell pair targets. Calibration rows excluded.',
      'objective':'Existing activity_loss: measured magnitude absolute error remains primary and unmodified by SE; same signed/RNA and paired auxiliaries. No detection labels, generated targets or AlphaGenome training features.',
      'optimization':{'head_lr':.0003,'reader_lr':.00001,'weight_decay':.001,'main_batch':256,'local_pairs':16,'cross_pairs':16,'gradient_clip':1.},
      'phase_one_gate':'Choose one common epoch by mean effect-size error skill across the three omitted neuronal sources. Adapted mean skill >=0.01, exceeds frozen control by >=0.005, no individual neuronal source error >1.05 times control, held-region nearby ordering and fixed-study/cell-pair rho each drop <=0.01, and nearby difference MAE beats zero-difference baseline. Otherwise stop this finite comparison.',
      'if_gate_passes':'Finish the other four whole-source roles at the same fixed three epochs, retaining the selected common epoch. Further model promotion needs all-source and outside-study evidence, useful confidence, and a separately frozen added-input test.',
      'stopping_rule':'Exactly this matched adaptation comparison; no epoch extension, learning-rate/loss-weight/architecture search, or new-test inspection to rescue failure.',
      'prior_work':'V8 already tried two epochs adapting the last convolutional stage on 150bp crops with older labels and an unidentified-fold checkpoint. This is not the first-ever core adaptation.',
      'new_test':'Ulirsch remains unused for model selection. Earlier blood/OL33/Matoba results remain revealed and immutable.',
      'scientific_goal_complete':False}
    path=DEST/'finite_native_adaptation_plan.json'
    if path.exists():assert json.loads(path.read_text())==plan
    else:write_json(path,plan)
    return plan
