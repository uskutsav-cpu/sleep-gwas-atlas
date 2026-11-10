"""Create a reviewable isolated trainer with one supervision change."""
from direction_common import *


def main():
    code=(MEASURE_WORK/'train_measurement_development.py').read_text()
    changes={
       'from measurement_common import *':'from direction_common import *',
       'from finite_measurement_plan import get_plan':'from finite_direction_plan import get_plan',
       'assert_prior_versions()':'assert_versions()',
       'd,pairs,bank,x,extra=load_data();train=':'d,pairs,bank,x,extra=load_data();new_signed=load_signed_sidecar(d);train=',
       "'new_Ulirsch_outcomes_used':False":"'new_original_signed_training_targets':int(new_signed[tr].sum()),'new_Ulirsch_outcomes_used':False",
       "sha(DEST/'finite_measurement_plan.json')":"sha(DEST/'finite_direction_plan.json')",
       "DEST/'folded_median_pretraining_audit.json'":"MEASUREMENT/'folded_median_pretraining_audit.json'",
       "Path(__file__).with_name('folded_effect_network.py')":"MEASURE_WORK/'folded_effect_network.py'",
       "Path(__file__).with_name('measurement_inputs.py')":"MEASURE_WORK/'measurement_inputs.py'",
    }
    for old,new in changes.items():
        assert old in code,old;code=code.replace(old,new)
    output=Path(__file__).with_name('train_direction_development.py')
    if output.exists():
        prior=json.loads((DEST/'trainer_derivation_receipt.json').read_text());assert sha(output)==prior['derived_source_sha256']
    output.write_text(code)
    write_json(DEST/'trainer_derivation_receipt.json',{'original_source_sha256':sha(MEASURE_WORK/'train_measurement_development.py'),
       'derived_source_sha256':sha(output),'explicit_replacements':changes,'only_training_data_change':'New source-verified signed labels; measured magnitudes, partitions, input features and actual batches unchanged.'})
    print('ISOLATED DIRECTION TRAINER WRITTEN')


if __name__=='__main__': main()
