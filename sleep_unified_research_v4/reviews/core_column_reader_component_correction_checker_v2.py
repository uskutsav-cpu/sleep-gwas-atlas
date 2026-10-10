"""Narrow initial-source/empty-input correction, only own synthetic fixtures."""
import ast
import hashlib
import importlib.util
import json
from pathlib import Path
import pickle
import sys
import pandas as pd

P=Path(__file__).resolve().parents[1];R=P/'reviews'
OLD=P/'scripts/core_column_reader_v1.py';NEW=P/'scripts/core_column_reader_v2.py'
ROOT=R/'core_column_reader_component_correction_controls_v2';ROOT.mkdir(exist_ok=False)
sha=lambda x:hashlib.sha256(Path(x).read_bytes()).hexdigest()
before={str(x):sha(x) for x in [OLD,NEW,R/'core_column_reader_component_review_seal_v1.json']}
assert before[str(NEW)]=='0aee12dc5f727be5a3242598a3cb4f7ba9ee11a7104ba257b821e4a2cf9c9cf4'
seal=json.loads((R/'core_column_reader_component_review_seal_v1.json').read_text())
for path,v in seal['artifacts'].items():assert sha(path)==v['sha256']
spec=importlib.util.spec_from_file_location('candidate_column_v2',NEW)
c=importlib.util.module_from_spec(spec);spec.loader.exec_module(c)
checks=[]
def check(name,condition,detail=None):
    checks.append({'control':name,'pass':bool(condition),'detail':detail})
    assert condition,name

# Exact whole-source transformation binds every method line, rather than
# treating version-name normalization as sufficient scientific equivalence.
expected=OLD.read_text().replace('def freeze_columns(path,read_kwargs,columns,spool,expected_rows):','def freeze_columns(path,read_kwargs,columns,spool,expected_rows,expected_source_sha256):')
expected=expected.replace("    python=read_kwargs.get('engine')=='python'","    if expected_rows<=0:raise RuntimeError('ORIGINAL_INPUT_HEADER_BUT_NO_ROWS_REJECTED')\n    initial_sha=digest(path)\n    if initial_sha!=expected_source_sha256:raise RuntimeError('SEALED_CORE_RAW_SOURCE_IDENTITY_DIFFERS_BEFORE_ANY_PARSER_READ')\n    python=read_kwargs.get('engine')=='python'")
expected=expected.replace("'source_sha256_before':digest(path),","'source_sha256_before':initial_sha,'expected_source_sha256':expected_source_sha256,")
expected=expected.replace("if result['source_sha256_before']!=result['source_sha256_after']:","if result['source_sha256_before']!=result['source_sha256_after'] or result['source_sha256_after']!=expected_source_sha256:")
check('exact_source_only_prespecified_correction',expected==NEW.read_text())
old_tree=ast.parse(OLD.read_text());new_tree=ast.parse(NEW.read_text())
for a,b in zip(old_tree.body,new_tree.body):
    if isinstance(a,ast.FunctionDef) and a.name!='freeze_columns':
        check('unchanged_method_AST_'+a.name,ast.dump(a,include_attributes=False)==ast.dump(b,include_attributes=False))

cols={'SNP':'snp','A1':'a1','A2':'a2','BETA':'beta'}
kw=lambda engine:{'sep':r'\s+' if engine=='python' else '\t','comment':None,'skiprows':0,'usecols':list(cols.values()),**({'engine':'python'} if engine=='python' else {'low_memory':True})}
initial='snp\ta1\ta2\tbeta\nrs1\ta\tc\t1\n';changed=initial.replace('\t1\n','\t9\n')

folder=ROOT/'former_false_PASS';folder.mkdir();path=folder/'source.txt';path.write_text(initial);wanted=sha(path)
gate=c.c_shape_gate
def mutate_after_shape(*args,**kwargs):
    result=gate(*args,**kwargs);path.write_text(changed);return result
c.c_shape_gate=mutate_after_shape
try:
    try:c.freeze_columns(path,kw('c'),cols,folder/'spool',1,wanted)
    except RuntimeError as e:check('former_shape_mutation_now_rejects','RAW_SOURCE_CHANGED' in str(e),str(e))
    else:raise AssertionError('Former witness still accepted')
finally:c.c_shape_gate=gate
check('former_witness_no_success_manifest',not (folder/'spool/column_preparation_manifest.json').exists())

for engine in ['c','python']:
    folder=ROOT/('wrong_initial_'+engine);folder.mkdir();path=folder/'source.txt'
    path.write_text(initial if engine=='c' else initial.replace('\t',' '))
    calls=[];gate=c.c_shape_gate;read=c.pd.read_csv
    def forbidden_shape(*a,**k):calls.append('shape');raise AssertionError('shape called before sealed identity')
    def forbidden_read(*a,**k):calls.append('parser');raise AssertionError('parser called before sealed identity')
    c.c_shape_gate=forbidden_shape;c.pd.read_csv=forbidden_read
    try:
        try:c.freeze_columns(path,kw(engine),cols,folder/'spool',1,'0'*64)
        except RuntimeError as e:check(engine+'_wrong_expected_reject','IDENTITY_DIFFERS_BEFORE_ANY_PARSER_READ' in str(e),str(e))
        else:raise AssertionError('Wrong expected source accepted')
    finally:c.c_shape_gate=gate;c.pd.read_csv=read
    check(engine+'_wrong_expected_no_source_shape_or_parser',calls==[])

for engine in ['c','python']:
    folder=ROOT/('empty_'+engine);folder.mkdir();path=folder/'source.txt'
    path.write_text('snp\ta1\ta2\tbeta\n' if engine=='c' else 'snp a1 a2 beta\n')
    calls=[];gate=c.c_shape_gate;read=c.pd.read_csv
    def forbidden_shape(*a,**k):calls.append('shape');raise AssertionError('shape called for globally empty source')
    def forbidden_read(*a,**k):calls.append('parser');raise AssertionError('parser called for globally empty source')
    c.c_shape_gate=forbidden_shape;c.pd.read_csv=forbidden_read
    try:
        try:c.freeze_columns(path,kw(engine),cols,folder/'spool',0,sha(path))
        except RuntimeError as e:check(engine+'_global_empty_reject','HEADER_BUT_NO_ROWS' in str(e),str(e))
        else:raise AssertionError('Empty original input accepted')
    finally:c.c_shape_gate=gate;c.pd.read_csv=read
    check(engine+'_empty_no_source_shape_or_parser',calls==[])

# Reuse two independent v1 accepted sources, binding every own fixture byte to
# the preserved map, then compare typed values exactly with the reviewed v1.
fixture_map=json.loads((R/'core_column_reader_component_fixture_identity_map_v1.json').read_text())['files']
for name in ['C_numeric_NA_bounds','Python_global_fallback']:
    oldfolder=R/'core_column_reader_component_controls_v1_3'/name
    path=oldfolder/'source.txt';old_manifest=json.loads((oldfolder/'spool/column_preparation_manifest.json').read_text())
    check(name+'_old_fixture_source_hash',sha(path)==fixture_map[str(path)]['sha256'])
    manifest=c.freeze_columns(path,old_manifest['read_kwargs'],old_manifest['columns'],ROOT/(name+'_spool'),old_manifest['expected_rows'],sha(path))
    check(name+'_before_after_expected_exact',manifest['source_sha256_before']==manifest['source_sha256_after']==manifest['expected_source_sha256']==sha(path))
    for target in old_manifest['typed_pieces']:
        def get_parts(m):
            out=[]
            for part in m['typed_pieces'][target]:
                assert sha(part['path'])==part['sha256']
                with open(part['path'],'rb') as f:out.append(pickle.load(f))
            return pd.concat(out)
        old=get_parts(old_manifest);new=get_parts(manifest)
        pd.testing.assert_series_equal(old,new,check_exact=True)
        check(name+'_'+target+'_unchanged_typed_value_dtype',True)
check('required_expected_source_argument',not ast.parse(NEW.read_text()).body[-1].args.defaults)
after={k:sha(k) for k in before};check('all_sources_and_prior_seal_unchanged',after==before)
for path,v in seal['artifacts'].items():assert sha(path)==v['sha256']
result={'status':'PASS_NARROW_COMPONENT_CORRECTION_PREPARATION_ONLY','checks':checks,'check_count':len(checks),'sources_before':before,'sources_after':after,'inherited_v1_check_count':468,'inherited_v1_seal_sha256':before[str(R/'core_column_reader_component_review_seal_v1.json')],'all18_inherited_review_artifacts_unchanged':True,'fixture_only':True,'execution_admitted':False,'full_adapter_supplied':False,'production_raw_body_reads':0,'reference_body_reads':0,'real_worker_or_fit_operations':0,'production_mutex_operations':0,'executable':sys.executable}
(R/'core_column_reader_component_correction_controls_receipt_v2.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'status':result['status'],'check_count':len(checks),'inherited_checks':468}))
