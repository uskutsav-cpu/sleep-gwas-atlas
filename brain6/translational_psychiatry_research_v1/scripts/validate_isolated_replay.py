"""Replay unchanged frozen calculations in scratch; preserve original hash-bound outputs.

Platform floating-point serialization can differ in final bits. Check every
replayed field numerically/string-wise, without overwriting the empirical record
or weakening any original source, code, QC or posterior check.
"""
from pathlib import Path
import hashlib,importlib.util,json,shutil,sys,tempfile,contextlib,io
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[3]
A=ROOT/'brain6/translational_psychiatry_research_v1'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
claims=pd.read_csv(A/'CLAIM_TO_EVIDENCE.tsv',sep='\t')
bindings={str(A/r.exact_input):sha(A/r.exact_input) for r in claims.itertuples()}
records=[]
with tempfile.TemporaryDirectory(prefix='brain6-replay-') as tmp:
    for name,folder,configs,tables in [
        ('global_contrasts','global_contrasts',['protocol.json'],['contrasts.tsv','covariance_sensitivity.tsv']),
        ('molecular_sensitivity','molecular',['sensitivity_protocol.json','serialization_amendment.json','coloc_5.2.3_source.R'],['sensitivity_results.tsv','independent_posterior_checks.tsv'])]:
        script=A/'scripts'/f'{name}.py';original=A/'research_v1'/folder;target=Path(tmp)/folder;target.mkdir()
        for file in configs:shutil.copy2(original/file,target/file)
        spec=importlib.util.spec_from_file_location(name,script);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
        module.OUT=target
        saved_args=sys.argv;sys.argv=[str(script),'run']
        try:
            with contextlib.redirect_stdout(io.StringIO()):module.main()
        finally:sys.argv=saved_args
        for table in tables:
            old=pd.read_csv(original/table,sep='\t');new=pd.read_csv(target/table,sep='\t')
            assert old.shape==new.shape and list(old.columns)==list(new.columns),table
            max_relative=0.;numeric=0
            for col in old:
                if pd.api.types.is_numeric_dtype(old[col]) and pd.api.types.is_numeric_dtype(new[col]):
                    a=old[col].to_numpy(dtype=float);b=new[col].to_numpy(dtype=float)
                    if col=='maximum_posterior_absolute_error':
                        # A floating-point discrepancy is judged by its original frozen
                        # replay tolerance, rather than relative equality to roundoff.
                        limit=json.loads((original/'sensitivity_protocol.json').read_text())['default_replay_tolerance']
                        assert np.isfinite(a).all() and np.isfinite(b).all()
                        assert (a>=0).all() and (b>=0).all() and max(a.max(),b.max())<=limit
                        numeric+=len(a)
                        continue
                    # Zero absolute tolerance retains relative precision for tiny P/BF values.
                    np.testing.assert_allclose(a,b,rtol=1e-12,atol=0,equal_nan=True,err_msg=table+':'+col)
                    mask=np.isfinite(a)&np.isfinite(b)&(a!=0)
                    if mask.any():max_relative=max(max_relative,float(np.max(abs((a[mask]-b[mask])/a[mask]))))
                    numeric+=len(a)
                else:assert old[col].fillna('').astype(str).equals(new[col].fillna('').astype(str)),table+':'+col
            records.append({'table':folder+'/'+table,'rows':len(old),'numeric_fields_checked':numeric,
                'maximum_relative_error':max_relative,'bitwise_equal':sha(original/table)==sha(target/table),'status':'PASS'})
after={path:sha(Path(path)) for path in bindings}
assert after==bindings,'Frozen result bytes were modified'
print(json.dumps({'classification':'UNCHANGED_CODE_SOURCE_VERIFIED_ISOLATED_NUMERICAL_REPLAY',
    'relative_tolerance':1e-12,'absolute_tolerance':0,'all_frozen_claim_artifacts_byte_unchanged':True,'tables':records},indent=2))
