"""Narrow synthetic-only v4 storage-change audit; no production source reads."""
import ast
import csv
import gc
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import pickle
import resource
import struct
import sys
import time
from types import SimpleNamespace
import warnings

import numpy as np
import pandas as pd

P=Path(__file__).resolve().parents[1]; R=P/'reviews'
PLAN=R/'independent_column_reader_v4_plan_v1.json'
plan=json.loads(PLAN.read_text()); root=Path(plan['synthetic_fixture_root'])
start=time.monotonic(); checks=[]; cases=[]; faults=[]
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(65536),b''):h.update(b)
    return h.hexdigest()
def check(name,condition,details=None):
    checks.append({'name':name,'pass':bool(condition),'details':details})
    if not condition:raise AssertionError(name)
def fence():
    limits=plan['bounds']
    if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss>limits['RSS_bytes']:raise RuntimeError('REVIEW_RSS_LIMIT')
    if time.monotonic()-start>limits['elapsed_seconds']:raise RuntimeError('REVIEW_DEADLINE')
    for path,key in [(P,'internal_emergency_bytes'),(root.parent,'SSD_floor_bytes')]:
        st=os.statvfs(path)
        if st.f_bavail*st.f_frsize<limits[key]:raise RuntimeError('REVIEW_FREE_SPACE_FLOOR')
    if root.exists() and sum(x.stat().st_size for x in root.rglob('*') if x.is_file())>limits['synthetic_namespace_bytes']:raise RuntimeError('REVIEW_NAMESPACE_CAP')
for path,value in plan['dependencies'].items():check('initial_identity_'+Path(path).name,sha(path)==value)
check('declared_runtime',sys.version_info[:3]==(3,11,11) and pd.__version__=='2.2.3' and np.__version__=='1.26.4')
fence();root.mkdir(exist_ok=False)
def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
v4=module('independent_candidate_v4',P/'scripts/core_column_reader_v4.py')
v3=module('independent_predecessor_v3',P/'scripts/core_column_reader_v3.py')
a=ast.parse((P/'scripts/core_column_reader_v3.py').read_text());b=ast.parse((P/'scripts/core_column_reader_v4.py').read_text())
old={n.name:n for n in a.body if isinstance(n,ast.FunctionDef)};new={n.name:n for n in b.body if isinstance(n,ast.FunctionDef)}
for name in old:
    if name!='freeze_columns':check('inherited_helper_AST_'+name,ast.dump(old[name])==ast.dump(new[name]))
check('inherited_nonfunction_module_AST',ast.dump(ast.Module(body=[n for n in a.body if not isinstance(n,ast.FunctionDef)],type_ignores=[]))==ast.dump(ast.Module(body=[n for n in b.body if not isinstance(n,ast.FunctionDef)],type_ignores=[])))
check('inherited_freeze_arguments',ast.dump(old['freeze_columns'].args)==ast.dump(new['freeze_columns'].args))
# Prefix through the unchanged store function, skipping only the docstring.
oldprefix=old['freeze_columns'].body[1:];newprefix=new['freeze_columns'].body[1:]
stop=next(i for i,n in enumerate(oldprefix) if isinstance(n,ast.FunctionDef) and n.name=='store')+1
check('inherited_initial_source_schema_shape_and_store_guards',all(ast.dump(x)==ast.dump(y) for x,y in zip(oldprefix[:stop],newprefix[:stop])) and len(newprefix[:stop])==stop)
check('inherited_50k_boundary',v4.ROWS==v3.ROWS==50000)
for path,value in json.loads((R/'core_column_reader_component_review_seal_v3.json').read_text())['artifacts'].items():
    check('sealed_v3_artifact_'+Path(path).name,sha(path)==value['sha256'] and Path(path).stat().st_size==value['bytes'])
main=next(n for n in ast.parse((P.parent/'scripts/01_harmonize.py').read_text()).body if isinstance(n,ast.FunctionDef) and n.name=='main')
nodes=[n for n in main.body if 417<=n.lineno<=454]
check('literal_original_AST_oracle_boundaries',nodes[0].lineno==417 and nodes[-1].lineno==454)
oracle=compile(ast.fix_missing_locations(ast.Module(body=nodes,type_ignores=[])),'literal_frozen_original_normalization','exec')
def fingerprint(s):
    h=hashlib.sha256(str(s.dtype).encode());h.update(np.asarray(s.isna(),dtype=np.uint8).tobytes())
    if isinstance(s.dtype,np.dtype) and s.dtype.kind in 'biufc':h.update(s.to_numpy(copy=False).tobytes())
    else:
        for value in s:
            h.update((type(value).__module__+'.'+type(value).__name__+':').encode())
            h.update(struct.pack('>d',float(value)) if isinstance(value,(float,np.floating)) else repr(value).encode());h.update(b'\0')
    return h.hexdigest()
def kwargs(engine,cols):
    d={'sep':'\t' if engine=='c' else r'\s+','comment':None,'skiprows':0,'usecols':list(cols.values())}
    d.update({'engine':'python'} if engine=='python' else {'low_memory':True});return d
def original(source,kw,cols,coord):
    raw=pd.read_csv(source,**kw)
    native={k:{'dtype':str(raw[k].dtype),'fingerprint':fingerprint(raw[k])} for k in raw}
    out=raw.rename(columns={v:k for k,v in cols.items()})
    env={'out':out,'pd':pd,'np':np,'log':lambda x:None,'fail':lambda x:(_ for _ in ()).throw(RuntimeError(x)),'args':SimpleNamespace(variant_map_strategy='BY_COORD_ALLELES' if coord else 'BY_RSID_ALLELES')}
    exec(oracle,env);return env['out'],native
def restore(m,target):
    parts=[]
    for info in m['typed_pieces'][target]:
        check('typed_hash_'+Path(info['path']).name,sha(info['path'])==info['sha256'])
        with Path(info['path']).open('rb') as f:s=pickle.load(f)
        check('typed_alignment_dtype_'+Path(info['path']).name,len(s)==info['rows'] and str(s.dtype)==info['dtype'] and s.index.equals(pd.RangeIndex(info['start'],info['start']+info['rows'])))
        parts.append(s)
    return pd.concat(parts)
def write_source(folder,cols,rows):
    folder.mkdir();path=folder/'synthetic.tsv'
    with path.open('x',newline='') as f:
        w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(cols.values());w.writerows(rows)
    return path
def run(name,cols,rows,engine='c',coord=False):
    if engine=='python':rows=([str(x).strip() for x in row] for row in rows)
    fence();folder=root/name;path=write_source(folder,cols,rows);wanted=sha(path);kw=kwargs(engine,cols)
    with warnings.catch_warnings(record=True) as ow:
        warnings.simplefilter('always');expected,native=original(path,kw,cols,coord)
    captures={};normal_calls=[];read=v4.pd.read_csv;tokens=v4.python_tokens;normalize=v4.normalized
    def capture_read(*args,**kw2):
        frame=read(*args,**kw2)
        if isinstance(frame,pd.DataFrame):
            for key in frame:captures[key]={'dtype':str(frame[key].dtype),'fingerprint':fingerprint(frame[key])}
        return frame
    def capture_tokens(*args,**kw2):
        for column,series,parts in tokens(*args,**kw2):
            captures[column]={'dtype':str(series.dtype),'fingerprint':fingerprint(series)}
            yield column,series,parts
    def capture_normalize(series,canonical):
        normal_calls.append({'canonical':canonical,'rows':len(series),'dtype':str(series.dtype),'index_start':int(series.index[0]) if len(series) else None})
        return normalize(series,canonical)
    v4.pd.read_csv=capture_read;v4.python_tokens=capture_tokens;v4.normalized=capture_normalize
    try:
        with warnings.catch_warnings(record=True) as nw:
            warnings.simplefilter('always');m4=v4.freeze_columns(path,kw,cols,folder/'v4_spool',len(expected),wanted,parse_coordinate_labels=coord)
    finally:v4.pd.read_csv=read;v4.python_tokens=tokens;v4.normalized=normalize
    check(name+'_full_native_inference_exact_dtypes_values_bits',captures==native,{'original':native,'candidate':captures})
    m3=v3.freeze_columns(path,kw,cols,folder/'v3_spool',len(expected),wanted,parse_coordinate_labels=coord)
    check(name+'_source_identity_before_after',m4['source_sha256_before']==m4['source_sha256_after']==sha(path)==wanted)
    check(name+'_target_dtypes_v3_equal',m4['global_dtype']==m3['global_dtype'])
    record={'name':name,'rows':len(expected),'engine':engine,'coordinate':coord,'source':str(path),'source_sha256':wanted,'native':native,'normalized':{},'original_warnings':[str(x.message) for x in ow],'candidate_warnings':[str(x.message) for x in nw]}
    for target in m4['typed_pieces']:
        got=restore(m4,target);prior=restore(m3,target)
        pd.testing.assert_series_equal(expected[target],got,check_names=False,check_exact=True)
        pd.testing.assert_series_equal(prior,got,check_names=False,check_exact=True)
        check(name+'_'+target+'_original_v3_v4_exact_dtype_NA_bits_index',fingerprint(expected[target])==fingerprint(prior)==fingerprint(got) and expected[target].index.equals(got.index))
        record['normalized'][target]={'dtype':str(got.dtype),'fingerprint':fingerprint(got),'missing':int(got.isna().sum()),'piece_rows':[x['rows'] for x in m4['typed_pieces'][target]]}
        if len(got)<20:record['normalized'][target]['values']=[repr(x) for x in got]
        del got,prior
    string_calls=[x for x in normal_calls if x['canonical'] in ['SNP','A1','A2']]
    check(name+'_elementwise_string_work_bounded',all(x['rows']<=50000 for x in string_calls))
    if coord:
        for target in ['CHR','BP']:
            check(name+'_'+target+'_numeric_inference_remains_global',len([x for x in normal_calls if x['canonical']==target])==1 and next(x for x in normal_calls if x['canonical']==target)['rows']==len(expected))
        check(name+'_coordinate_tokens_spooled_string_indices',all(x['dtype']=='string' and x['rows']<=50000 for parts in m4['coordinate_token_pieces'].values() for x in parts))
    cases.append(record);del expected,m3,m4;gc.collect();fence()

small={'SNP':'label','A1':'allele1','A2':'allele2','BETA':'beta','N':'n'}
sets=[('integer_native',[['1','1','0','-0.0','001'],['9223372036854775807','0','1','1.00000000000000001','2']]),
      ('float_native',[['-0.0','1e-7','1.2345678901234567','-0.0','NA'],['1e20','NaN','Inf','1.0000000000000002','3']]),
      ('bool_native',[['True','False','TRUE','1','4'],['False','True','FALSE','2','5']]),
      ('mixed_object',[['True','False','a','-1','NA'],['001','1e0','NULL','.1','True'],['rsMixed','  t  ','g','bad','False']]),
      ('all_missing',[['NA','NULL','NaN','NaN','NA'],['NULL','NaN','NA','NA','NULL']])]
for engine in ['c','python']:
    for name,rows in sets:run(engine+'_'+name,small,rows,engine)

coord={'SNP':'label','A1':'allele1','A2':'allele2','BETA':'beta'}
labels=[' chr01:00042:A:T ','CHR1:42:A:T','Chr2:3','2_0003_G_A','chrX:4','chr1.0:2','chr1:2.0','chr1:2tail','NA','001','1e0','chr18446744073709551616:42','chr2:18446744073709551616','chr3:'+('9'*310),'chr1:9223372036854775807','chr2:9223372036854775808','chr3:18446744073709551615']
for engine in ['c','python']:
    run(engine+'_coordinate_extremes_case_missing',coord,[[x,'a','c','-0.0'] for x in labels],engine,True)
    run(engine+'_coordinate_all_NA',coord,[['NA','NA','NULL','1'],['CHR1:2','a','c','1']],engine,True)
    # A decisive late value changes the dtype of the complete source column.
    rows=(['001' if i<50000 else ['True','rsMiXeD','NA'][i-50000], 'True' if i<50000 else ['001','a','NA'][i-50000], '1.0000000000000002' if i<50000 else ['False','t','NA'][i-50000], '-0.0' if i%2==0 else '1.0000000000000002','001'] for i in range(50003))
    run(engine+'_strings_full_inference_cross_50k',small,rows,engine,False)
    rows=(['chr01:00042:A:T' if i<50000 else ['CHR1:42','chr2:18446744073709551616','NA'][i-50000], 'a','c','-0.0'] for i in range(50003))
    run(engine+'_coordinate_full_nullable_overflow_cross_50k',coord,rows,engine,True)

def reject(name,engine,hook=None,wanted_override=None):
    fence();folder=root/name;path=write_source(folder,coord,[['chr1:42','a','c','1'],['chr2:43','t','g','2']]);wanted=sha(path)
    saved={key:getattr(v4,key) for key in ['normalized','dump_new','digest','load']}
    if hook:hook(folder,path,saved)
    try:
        try:v4.freeze_columns(path,kwargs(engine,coord),coord,folder/'v4_spool',2,wanted if wanted_override is None else wanted_override,True)
        except RuntimeError as e:error=str(e)
        else:raise AssertionError(name+'_unexpected_accept')
    finally:
        for key,value in saved.items():setattr(v4,key,value)
    check(name+'_reject_no_success_manifest',not (folder/'v4_spool/column_preparation_manifest.json').exists(),error)
    faults.append({'name':name,'error':error});return error
for engine in ['c','python']:
    check(engine+'_wrong_initial_sha_fail_closed','IDENTITY_DIFFERS_BEFORE_ANY_PARSER_READ' in reject(engine+'_initial_source_hash',engine,wanted_override='0'*64))
    def change_source(folder,path,saved):
        once=[]
        def normalization(series,target):
            if not once:path.write_text(path.read_text().replace('chr1:42','chr1:99'));once.append(True)
            return saved['normalized'](series,target)
        v4.normalized=normalization
    check(engine+'_midwork_source_sha_fail_closed','RAW_SOURCE_CHANGED' in reject(engine+'_midwork_source_change',engine,change_source))
    def change_token(folder,path,saved):
        once=[]
        def normalization(series,target):
            if target=='SNP' and not once:
                token=next((folder/'v4_spool').glob('coordinate_tokens_CHR_*.pkl'))
                with token.open('ab') as f:f.write(b'changed')
                once.append(True)
            return saved['normalized'](series,target)
        v4.normalized=normalization
    check(engine+'_coordinate_token_hash_fail_closed','COORDINATE_TOKEN_SPOOL_CHANGED' in reject(engine+'_token_hash_change',engine,change_token))
    for variant in ['index','dtype','length']:
        def invalid_token(folder,path,saved,variant=variant):
            def dump(dest,value):
                if Path(dest).name.startswith('coordinate_tokens_CHR_'):
                    value=value.copy()
                    if variant=='index':value.index=value.index+1
                    elif variant=='dtype':value=value.astype(object)
                    else:value=value.iloc[:-1]
                return saved['dump_new'](dest,value)
            v4.dump_new=dump
        check(engine+'_'+variant+'_token_alignment_fail_closed','COORDINATE_TOKEN_ALIGNMENT_OR_DTYPE_CHANGED' in reject(engine+'_token_'+variant,engine,invalid_token))
fence()
for path,value in plan['dependencies'].items():check('final_identity_'+Path(path).name,sha(path)==value)
files={str(x):{'bytes':x.stat().st_size,'sha256':sha(x)} for x in sorted(root.rglob('*')) if x.is_file()}
with (R/'independent_column_reader_v4_fixture_map_v1.json').open('x') as f:json.dump({'invented_fixtures_only':True,'files':files},f,indent=2);f.write('\n')
receipt={'status':'PASS_COMPONENT_PREPARATION_ONLY_NOT_EXECUTION_ADMISSION','checks':checks,'check_count':len(checks),'cases':cases,'case_count':len(cases),'faults':faults,'fault_count':len(faults),'plan_sha256':sha(PLAN),'dependencies_before_after':plan['dependencies'],'prior_v3_517_controls_inherited_not_rerun':True,'elapsed_seconds':time.monotonic()-start,'peak_RSS_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'synthetic_namespace_bytes':sum(v['bytes'] for v in files.values()),'runtime':{'python':sys.version,'numpy':np.__version__,'pandas':pd.__version__,'executable':sys.executable},'production_GWAS_reference_reads':0,'production_mutex_operations':0,'production_worker_replay_estimator_calls':0,'actual_2GiB_full_source_memory_certified':False,'execution_admitted':False}
with (R/'independent_column_reader_v4_receipt_v1.json').open('x') as f:json.dump(receipt,f,indent=2);f.write('\n')
print(json.dumps({key:receipt[key] for key in ['status','check_count','case_count','fault_count','elapsed_seconds','peak_RSS_bytes','synthetic_namespace_bytes']}))
