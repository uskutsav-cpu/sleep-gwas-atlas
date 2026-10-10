"""Small actual-frame diagnosis only; no raw body or production worker."""
import ast, hashlib, importlib.util, json, pathlib, pickle, sys
import pandas as pd
import numpy as np
P=pathlib.Path(__file__).resolve().parents[1];sys.path.insert(0,str(P/'scripts'))
import core_bounded_harmonizer_v3 as old
import core_bounded_harmonizer_v4 as new
original=P.parent/'scripts/01_harmonize.py'
module,body=old.load_original(original)
oldst=old.original_stages(module,body);newst=new.original_stages(module,body)
folder=pathlib.Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4/core_pipeline/large35_checkpoint_continuation_v1/sleep_apnea/ephemeral_bounded_spool/post_liftover')
observed={};cases=[]
for index in [388,389]:
    path=folder/('%08d.pkl'%index);payload=path.read_bytes();frame=pickle.loads(payload)
    observed[str(path)]={'sha256':hashlib.sha256(payload).hexdigest(),'bytes':len(payload),'rows':len(frame),'dtypes':{k:str(v) for k,v in frame.dtypes.items()}}
    if index==389:
        assert len(frame)==0
        try:oldst['qc'](frame.copy(),lambda s:s.duplicated(keep='first'))
        except TypeError as error:
            assert str(error)=="Cannot perform reduction 'sum' with string dtype"
            cases.append({'case':'ACTUAL_EMPTY389_ORIGINAL_REDUCTION_FAILURE','status':'EXPECTED_FAILURE_PRESERVED','reason':str(error)})
        else:raise AssertionError('Actual old failure was not reproduced')
        out,steps=newst['qc'](frame.copy(),lambda s:s.duplicated(keep='first'))
        assert len(out)==0 and all(removed==remaining==0 for _,removed,remaining in steps)
        assert len(steps)==13
        observed['corrected_empty_ordinary_steps']=steps
        cases.append({'case':'ZERO_ROW_ONLY_DROP_REPAIR_ALL13_ORDERED_ORDINARY_STEPS','status':'PASS'})
    else:
        assert len(frame)==42291
        before,steps_before=oldst['qc'](frame.copy(),lambda s:s.duplicated(keep='first'))
        after,steps_after=newst['qc'](frame.copy(),lambda s:s.duplicated(keep='first'))
        pd.testing.assert_frame_equal(before,after,check_exact=True)
        assert steps_before==steps_after
        observed['adjacent_nonempty_ordinary_steps']=steps_before
        cases.append({'case':'ACTUAL_NONEMPTY388_ORIGINAL_AND_CORRECTED_ROWS_DTYPES_BITS_QC_EQUAL','status':'PASS'})
        # Exact same dataframe, first predicate retained. A later filter empties
        # the frame, so the subsequent rsID/allele/missingness reductions must
        # also remain zero. This is a representative actual-frame slice.
        later=frame.loc[frame.index[:3]].copy();later['A1']=later['A2']
        out,steps=newst['qc'](later,lambda s:s.duplicated(keep='first'))
        assert len(out)==0
        at=next(i for i,x in enumerate(steps) if x[0]=='A1 equals A2')
        assert all(removed==remaining==0 for _,removed,remaining in steps[at+1:])
        cases.append({'case':'ACTUAL_ROW_SLICE_LATER_EMPTY_MASKS_ACCOUNT_ZERO','status':'PASS'})
# Boundaries remain original; only stage_qc drop reducer is amended.
for name in ['life','required','mapping','sample','serialize']:
    assert oldst[name].__code__.co_code==newst[name].__code__.co_code
    assert oldst[name].__code__.co_consts==newst[name].__code__.co_consts
for name in ['load_original','frames','aggregate_steps','save_frame','compiled_stage']:
    a=next(n for n in ast.parse((P/'scripts/core_bounded_harmonizer_v3.py').read_text()).body if isinstance(n,ast.FunctionDef) and n.name==name)
    b=next(n for n in ast.parse((P/'scripts/core_bounded_harmonizer_v4.py').read_text()).body if isinstance(n,ast.FunctionDef) and n.name==name)
    assert ast.dump(a,include_attributes=False)==ast.dump(b,include_attributes=False)
receipt={'status':'ACTUAL_EMPTY_CHUNK_CAUSE_CONFIRMED_ZERO_ROW_CORRECTION_FOCUSED_PASS','runtime':{'pandas':pd.__version__,'numpy':np.__version__},'cases':cases,'case_count':len(cases),'actual_small_checkpoint_evidence':observed,'original_code_sha256':hashlib.sha256(original.read_bytes()).hexdigest(),'old_helper_sha256':hashlib.sha256((P/'scripts/core_bounded_harmonizer_v3.py').read_bytes()).hexdigest(),'corrected_helper_sha256':hashlib.sha256((P/'scripts/core_bounded_harmonizer_v4.py').read_bytes()).hexdigest(),'source_raw_body_reads':0,'source_reparsed':False,'workers_launched':0,'production_mutex_operations':0,'old_suite_repeated':False}
path=P/'reviews/core_apnea_empty_frame_diagnosis_controls_receipt_v1.json'
with path.open('x') as f:json.dump(receipt,f,indent=2);f.write('\n')
print(json.dumps({'receipt':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'cases':cases},indent=2))
