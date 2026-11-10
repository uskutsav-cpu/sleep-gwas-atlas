import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'sleep_frailty_ml_direction'))
from direction_common import *
DIRECTION=DEST
DEST=ROOT/'outputs/sleep_frailty_ml_representation'


def assert_versions_all():
    assert_versions()
    for name,h in json.loads((DIRECTION/'final_progress_verification.json').read_text())['files'].items():
        assert sha(DIRECTION/name)==h,name


def get_plan():
    plan={'purpose':'Determine whether fixed DNA features contain a useful simple effect-size signal before another neural architecture experiment.',
       'data_sha256':sha(OLD/'coverage_training_arrays.npz'),'signed_sidecar_sha256':sha(DIRECTION/'original_paired_signed_development.npy'),
       'roles':['none','Salomon','Kosicki2025','Siraj','IGVF240K','Koesterich2023','McAfee2023','Retallick2026'],
       'magnitude_features':'96 absolute projected channel differences; 42 symmetric scalar reader summaries; actual length and changed position. Add six scalar readout-by-supported-line interactions and line offsets.',
       'signed_features':'96 signed channel differences plus three signed readout differences; three readout-by-supported-line interactions. No intercept or centered differences, preserving allele-swap antisymmetry.',
       'normalization':'Unique physical training banks only; training mean and SD for magnitude features, RMS without centering for signed features. Clip normalized features at 12.',
       'magnitude_objective':'Weighted mean sqrt((predicted actual magnitude - observed magnitude)^2 + 0.01^2), plus 0.1/2 times squared standardized slope norm. Line offsets are unpenalized. No SE weighting of primary magnitude error.',
       'signed_objective':'Weighted mean squared original signed-effect residual plus 0.1 times squared slope norm; documented source signed labels only; same training-derived context and precision weights as previous signed auxiliary.',
       'optimization':'Magnitude: one L-BFGS-B fit per role, maximum 80 iterations and 120 function evaluations. Signed: one closed-form ridge solve per role. Fixed penalties; no grid or model/epoch search.',
       'endpoint':'Actual magnitude MAE versus training-only typical line/cell effect, macro across supported cells within each held source. Report nearby and fixed-study/cell-pair differences separately.',
       'signed_output_scope':'Separate diagnostic of direction learnability; abs(signed prediction) is not substituted for the primary observed magnitude prediction.',
       'safety':'Clip negative magnitude predictions to zero; no effect-size, significance or SE-based test-row filtering. No new outside study, calibrated confidence, project ranks or patient-risk claim.',
       'stopping_rule':'Finish exactly these two simple readout diagnostics. A negative result does not prove the representation contains no nonlinear signal. Do not tune more penalties or open Ulirsch to find a favorable result.',
       'AlphaGenome_training_inputs':False,'pretrained_weights_updated':False,'reserved_study_used':False,'scientific_success_assumed':False}
    path=DEST/'finite_representation_plan.json'
    if path.exists(): assert json.loads(path.read_text())==plan
    else: write_json(path,plan)
    return plan
