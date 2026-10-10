"""Independent component fixtures only: no production inputs or execution."""
import ast
import csv
import gzip
import hashlib
import importlib.util
import json
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

P=Path(__file__).resolve().parents[1];R=P/'reviews'
NEW=P/'scripts/core_column_reader_v3.py';OLD=P/'scripts/core_column_reader_v2.py'
ORIGINAL=P.parent/'scripts/01_harmonize.py'
ROOT=R/'core_column_reader_component_controls_v3_1';ROOT.mkdir(exist_ok=False)
CHECKS=[];CASES=[];START=time.monotonic()
sha=lambda path:hashlib.sha256(Path(path).read_bytes()).hexdigest()
def check(name,condition,detail=None):
    CHECKS.append({'control':name,'pass':bool(condition),'detail':detail})
    if not condition:raise AssertionError(name)

prior=json.loads((R/'core_column_reader_component_correction_review_v2.json').read_text())
BINDS=[NEW,OLD,ORIGINAL,R/'core_column_reader_component_review_seal_v1.json',R/'core_column_reader_component_correction_review_seal_v2.json']+[Path(x) for x in prior['inherited_runtime_code_identities_still_match']]
BEFORE={str(x):sha(x) for x in BINDS}
check('candidate_exact_identity',BEFORE[str(NEW)]=='6c0ecebce59af3f73986e059f0facd0329b9168f9fb89e2e78db39ca543709b0')
check('original_exact_identity',BEFORE[str(ORIGINAL)]=='487074007f8361e7de7f8f771dd0a912f926d4b7a0427cfa3c47b02be1421f7a')
for path,value in prior['inherited_runtime_code_identities_still_match'].items():check('inherited_identity_'+path,sha(path)==value)
for name in ['core_column_reader_component_review_seal_v1.json','core_column_reader_component_correction_review_seal_v2.json']:
    seal=json.loads((R/name).read_text())
    for path,value in seal['artifacts'].items():check('prior_sealed_artifact_'+Path(path).name,sha(path)==value['sha256'] and Path(path).stat().st_size==value['bytes'])
check('declared_runtime',sys.version_info[:3]==(3,11,11) and np.__version__=='1.26.4' and pd.__version__=='2.2.3')
spec=importlib.util.spec_from_file_location('candidate_v3',NEW);c=importlib.util.module_from_spec(spec);spec.loader.exec_module(c)
oldtree=ast.parse(OLD.read_text());newtree=ast.parse(NEW.read_text())
oldfunc={x.name:x for x in oldtree.body if isinstance(x,ast.FunctionDef)}
newfunc={x.name:x for x in newtree.body if isinstance(x,ast.FunctionDef)}
for name in oldfunc:
    if name!='freeze_columns':check('unchanged_helper_AST_'+name,ast.dump(oldfunc[name],include_attributes=False)==ast.dump(newfunc[name],include_attributes=False))
oldfreeze=oldfunc['freeze_columns'];newfreeze=newfunc['freeze_columns']
check('unchanged_initial_guards_before_and_through_shape',all(ast.dump(a,include_attributes=False)==ast.dump(b,include_attributes=False) for a,b in zip(oldfreeze.body[1:12],newfreeze.body[1:12])))
check('optional_coordinate_default_false',len(newfreeze.args.defaults)==1 and isinstance(newfreeze.args.defaults[0],ast.Constant) and newfreeze.args.defaults[0].value is False)
check('expected_source_still_required',newfreeze.args.args[-2].arg=='expected_source_sha256')
check('spool_block_size_unchanged',c.ROWS==50000)

# Execute the exact original branch and whole-normalization AST. This oracle
# starts from one full original selected-column read, never the candidate.
tree=ast.parse(ORIGINAL.read_text());main=next(x for x in tree.body if isinstance(x,ast.FunctionDef) and x.name=='main')
nodes=[x for x in main.body if 416<=x.lineno<=454]
check('literal_original_oracle_range',nodes[0].lineno==417 and nodes[-1].lineno==454)
ORACLE=compile(ast.fix_missing_locations(ast.Module(body=nodes,type_ignores=[])),str(ORIGINAL),'exec')
def fail(message):raise RuntimeError(message)
def fingerprint(s):
    h=hashlib.sha256(str(s.dtype).encode());h.update(np.asarray(s.isna(),dtype=np.uint8).tobytes())
    if isinstance(s.dtype,np.dtype) and s.dtype.kind in 'biufc':h.update(s.to_numpy(copy=False).tobytes())
    else:
        for value in s:
            h.update((type(value).__module__+'.'+type(value).__name__+':').encode())
            h.update(struct.pack('>d',float(value)) if isinstance(value,(float,np.floating)) else repr(value).encode());h.update(b'\0')
    return h.hexdigest()
def kw(engine,columns,sep=None,skip=0):
    result={'sep':sep or (r'\s+' if engine=='python' else '\t'),'comment':None,'skiprows':skip,'usecols':list(dict.fromkeys(columns.values()))}
    result.update({'engine':'python'} if engine=='python' else {'low_memory':True});return result
def original(path,kwargs,columns,coord):
    raw=pd.read_csv(path,**kwargs)
    parsed={x:{'dtype':str(raw[x].dtype),'fingerprint':fingerprint(raw[x])} for x in raw}
    builds=[str(x) for x in raw[columns['BUILD']].dropna().unique() if str(x).strip()] if 'BUILD' in columns else []
    out=raw.rename(columns={v:k for k,v in columns.items()})
    env={'out':out,'pd':pd,'np':np,'log':lambda x:None,'fail':fail,'args':SimpleNamespace(variant_map_strategy='BY_COORD_ALLELES' if coord else 'BY_RSID_ALLELES')}
    exec(ORACLE,env);return env['out'],parsed,builds,env['coordinate_label_parsed']
def restore(manifest,target):
    out=[]
    for part in manifest['typed_pieces'][target]:
        check(target+'_typed_hash_'+Path(part['path']).name,sha(part['path'])==part['sha256'])
        with open(part['path'],'rb') as f:s=pickle.load(f)
        check(target+'_typed_dtype_index_'+Path(part['path']).name,str(s.dtype)==part['dtype'] and s.index.equals(pd.RangeIndex(part['start'],part['start']+part['rows'])))
        out.append(s)
    return pd.concat(out)
def run(name,text,columns,engine='c',coord=True,sep=None,compressed=False,skip=0):
    folder=ROOT/name;folder.mkdir();path=folder/('source.txt.gz' if compressed else 'source.txt')
    if compressed:
        with gzip.open(path,'wb') as f:f.write(text.encode())
    else:path.write_text(text)
    kwargs=kw(engine,columns,sep,skip)
    with warnings.catch_warnings(record=True) as oldwarnings:
        warnings.simplefilter('always');out,parsed,builds,original_flag=original(path,kwargs,columns,coord)
    inferred={};normalizer=c.normalized
    def capture(series,canonical):
        if canonical in columns:inferred[columns[canonical]]={'dtype':str(series.dtype),'fingerprint':fingerprint(series)}
        return normalizer(series,canonical)
    c.normalized=capture
    try:
        with warnings.catch_warnings(record=True) as newwarnings:
            warnings.simplefilter('always');manifest=c.freeze_columns(path,kwargs,columns,folder/'spool',len(out),sha(path),parse_coordinate_labels=coord)
    finally:c.normalized=normalizer
    check(name+'_full_source_parser_exact_types_and_bits',inferred==parsed,{'original':parsed,'candidate':inferred})
    check(name+'_source_initial_final_expected',manifest['source_sha256_before']==manifest['source_sha256_after']==manifest['expected_source_sha256']==sha(path))
    check(name+'_global_BUILD_values',manifest['original_BUILD_values']==builds)
    check(name+'_coordinate_flag',manifest['coordinate_label_parsed']==original_flag==coord)
    targets={('BETA' if x=='OR' else x) for x in columns}|({'CHR','BP'} if coord else set())
    check(name+'_exact_target_set',set(manifest['typed_pieces'])==targets)
    record={'name':name,'rows':len(out),'source_sha256':sha(path),'source_bytes':path.stat().st_size,'kwargs':kwargs,'columns':columns,'parse_coordinate_labels':coord,'normalized':{},'original_warnings':[str(x.message) for x in oldwarnings],'candidate_warnings':[str(x.message) for x in newwarnings]}
    for target in sorted(targets):
        restored=restore(manifest,target);pd.testing.assert_series_equal(out[target],restored,check_names=False,check_exact=True)
        check(name+'_'+target+'_global_dtype_value_NA_IEEE_bits',str(out[target].dtype)==manifest['global_dtype'][target] and fingerprint(out[target])==fingerprint(restored))
        record['normalized'][target]={'dtype':str(restored.dtype),'fingerprint':fingerprint(restored),'NA_count':int(restored.isna().sum()),'parts':len(manifest['typed_pieces'][target])}
        if len(out)<=25:record['normalized'][target]['literal_values']=[repr(x) for x in restored]
    if engine=='c':
        widths=[];blocks=[]
        for selected in [kwargs['usecols'],[kwargs['usecols'][0]]]:
            reader=pd.read_csv(path,**dict(kwargs,usecols=selected,iterator=True));widths.append(int(reader._engine._reader.table_width))
            values=reader._engine._reader.read_low_memory(None);blocks.append([len(next(iter(x.values()))) for x in values]);reader.close()
        check(name+'_C_physical_width_and_native_blocks',widths[0]==widths[1]==manifest['physical_shape_gate']['full_header_width'] and blocks[0]==blocks[1])
        record['native_C_widths']=widths;record['native_C_blocks']=blocks
    CASES.append(record);return manifest,out

BASE={'SNP':'label','A1':'a1','A2':'a2','BETA':'beta','SE':'se','P':'p','N':'n','BUILD':'build'}
labels=['chr1:00042:A:T','CHR1:00042:A:T','Chr2:3','2_0003_G_A','chr0003:00004','4:5','chr5_6_x','chr1.0:2','chr1:2.0','chr6:7tail','NA','001','1e0','7:8_','chr8:9:extra','chrX:12','chr9:10:-']
header=list(BASE.values())
rows=[[label,'a','c','-0.0' if i==0 else ('bad' if i==8 else '1.00000000000000001'),'1e-2','.5','True' if i==0 else ('False' if i==1 else '001'),'hg19' if i%2 else 'GRCh37'] for i,label in enumerate(labels)]
for engine in ['c','python']:
    sep='\t' if engine=='c' else ' '
    text=sep.join(header)+'\n'+''.join(sep.join(row)+'\n' for row in rows)
    m,out=run(engine+'_coordinate_literal_cases',text,BASE,engine=engine)
    check(engine+'_uppercase_prefix_still_invalid_after_SNP_lower',out['SNP'].iloc[1]=='chr1:00042:a:t' and pd.isna(out['CHR'].iloc[1]) and out['CHR'].iloc[0]==1 and out['BP'].iloc[0]==42)
    check(engine+'_mixed_labels_full_nullable_Int64',str(out['CHR'].dtype)=='Int64' and str(out['BP'].dtype)=='Int64')
    check(engine+'_decimal_and_suffix_boundaries_invalid',out['CHR'].iloc[[7,8,9]].isna().all())

MIN={'SNP':'label','A1':'a1','A2':'a2','OR':'odds'}
run('C_comma_stripped_labels_OR_gzip','metadata\nlabel,a1,a2,odds,unused\n" chr01:0002:A:T "," a ",c,1.0000000000000002,"a,b"\nCHR2:3,t,g,0,x\nNA,a,c,-1,x\n',MIN,sep=',',skip=1,compressed=True)
for engine in ['c','python']:
    sep='\t' if engine=='c' else ' '
    run(engine+'_all_missing_labels',sep.join(['label','a1','a2','odds'])+'\n'+sep.join(['NA','NA','NULL','NaN'])+'\n'+sep.join(['NULL','a','c','0'])+'\n',MIN,engine=engine)
    bounds=['chr1:9223372036854775807','chr2:9223372036854775808','chr3:18446744073709551615','chr4:18446744073709551616','NA','chr18446744073709551616:42']
    run(engine+'_coordinate_uint64_overflow',sep.join(['label','a1','a2','odds'])+'\n'+''.join(sep.join([label,'a','c','1'])+'\n' for label in bounds),MIN,engine=engine)
run('Python_short_blank_rows','label a1 a2 odds\nchr1:1 a c 1\n\nNA\nchr2:2 t g\nchr3:3 a\n',MIN,engine='python')

# Physical width changes C low_memory block boundaries. Both full original
# parsing and singleton candidate parsing see width64; typed pieces cross50k.
widehead=['label','a1','a2','beta','n','build']+['unused'+str(i) for i in range(58)]
for engine in ['c','python']:
    sep='\t' if engine=='c' else ' ';text=sep.join(widehead)+'\n'
    lines=[]
    for i in range(50003):
        label='chr1:00042:A:T' if i<50000 else ['CHR1:42','chr1:42.0','2_00043'][i-50000]
        beta='True' if i<8192 else ('001' if i<16384 else ('1.00000000000000001' if i<50000 else 'bad'))
        lines.append(sep.join([label,'a','c',beta,'001','hg19' if i<50000 else 'GRCh37']+['0']*58)+'\n')
    manifest,out=run(engine+'_coordinate_cross_50k_and_full_width',text+''.join(lines),{'SNP':'label','A1':'a1','A2':'a2','BETA':'beta','N':'n','BUILD':'build'},engine=engine)
    check(engine+'_two_typed_chunks_global_coordinate_dtype',str(out['CHR'].dtype)=='Int64' and str(out['BP'].dtype)=='Int64' and [x['rows'] for x in manifest['typed_pieces']['CHR']]==[50000,3])
    del lines,text,manifest,out

# Retention is a storage-lifetime change. Reuse four previously sealed own
# accepted sources and compare v3 false/default flag to the literal original.
oldmap=json.loads((R/'core_column_reader_component_fixture_identity_map_v1.json').read_text())['files']
for name in ['C_numeric_NA_bounds','C_OR_full_log','Python_global_fallback','Python_bool_numeric']:
    folder=R/'core_column_reader_component_controls_v1_3'/name
    paths=list(folder.glob('source.txt*'));source=paths[0];m=json.loads((folder/'spool/column_preparation_manifest.json').read_text())
    check(name+'_sealed_fixture_reuse',sha(source)==oldmap[str(source)]['sha256'])
    text=gzip.decompress(source.read_bytes()).decode() if source.name.endswith('.gz') else source.read_text()
    run('retention_'+name,text,m['columns'],engine='python' if m['read_kwargs'].get('engine')=='python' else 'c',coord=False,sep=m['read_kwargs']['sep'],compressed=source.name.endswith('.gz'))

def reject(name,text,columns,engine='c',coord=True,rows=1,expected=None,patch=None,error=''):
    folder=ROOT/name;folder.mkdir();path=folder/'source.txt';path.write_text(text);kwargs=kw(engine,columns)
    if patch:patch(path)
    try:c.freeze_columns(path,kwargs,columns,folder/'spool',rows,sha(path) if expected is None else expected,parse_coordinate_labels=coord)
    except RuntimeError as e:check(name+'_reject',error in str(e),str(e))
    else:raise AssertionError(name+'_accepted')
    check(name+'_no_success_manifest',not (folder/'spool/column_preparation_manifest.json').exists())
for engine in ['c','python']:
    sep='\t' if engine=='c' else ' '
    for suffix,cols in [('no_SNP',{'A1':'a1','A2':'a2','BETA':'beta'}),('only_CHR',{'SNP':'label','CHR':'chrom','A1':'a1','A2':'a2','BETA':'beta'}),('only_BP',{'SNP':'label','BP':'pos','A1':'a1','A2':'a2','BETA':'beta'}),('both_CHR_BP',{'SNP':'label','CHR':'chrom','BP':'pos','A1':'a1','A2':'a2','BETA':'beta'})]:
        reject(engine+'_coordinate_branch_guard_'+suffix,sep.join(cols.values())+'\n'+sep.join(['1']*len(cols))+'\n',cols,engine,error='EXACT_ORIGINAL_COORDINATE_LABEL_BRANCH_REQUIRED')
    cols={'SNP':'label','A1':'a1','A2':'a2','BETA':'beta'}
    text=sep.join(cols.values())+'\n'+sep.join(['chr1:1','a','c','1'])+'\n'
    calls=[];read=c.pd.read_csv;shape=c.c_shape_gate
    def forbid(*a,**k):calls.append('parser/shape');raise AssertionError('identity/empty gate did not precede parsing')
    c.pd.read_csv=forbid;c.c_shape_gate=forbid
    try:
        reject(engine+'_wrong_sealed_identity',text,cols,engine,expected='0'*64,error='IDENTITY_DIFFERS_BEFORE_ANY_PARSER_READ')
        reject(engine+'_empty_reject',sep.join(cols.values())+'\n',cols,engine,rows=0,error='HEADER_BUT_NO_ROWS')
    finally:c.pd.read_csv=read;c.c_shape_gate=shape
    check(engine+'_identity_and_empty_preparser',calls==[])
folder=ROOT/'C_postshape_source_mutation';folder.mkdir();path=folder/'source.txt';path.write_text('label\ta1\ta2\tbeta\nchr1:1\ta\tc\t1\n');wanted=sha(path);shape=c.c_shape_gate
def mutate(*args,**kwargs):
    result=shape(*args,**kwargs);path.write_text(path.read_text().replace('\t1\n','\t9\n'));return result
c.c_shape_gate=mutate
try:
    try:c.freeze_columns(path,kw('c',{'SNP':'label','A1':'a1','A2':'a2','BETA':'beta'}),{'SNP':'label','A1':'a1','A2':'a2','BETA':'beta'},folder/'spool',1,wanted,parse_coordinate_labels=True)
    except RuntimeError as e:check('postshape_identity_mutation_still_rejects','RAW_SOURCE_CHANGED' in str(e),str(e))
    else:raise AssertionError('postshape_source_mutation_accepted')
finally:c.c_shape_gate=shape
check('postshape_mutation_no_success_manifest',not (folder/'spool/column_preparation_manifest.json').exists())
AFTER={path:sha(path) for path in BEFORE};check('all_bound_sources_before_after_unchanged',BEFORE==AFTER)
receipt={'status':'PASS_V3_COMPONENT_PREPARATION_ONLY','checks':CHECKS,'check_count':len(CHECKS),'cases':CASES,'case_count':len(CASES),'sources_before':BEFORE,'sources_after':AFTER,'inherited_v1_controls':468,'inherited_v2_controls':40,'executable':sys.executable,'runtime':{'python':sys.version,'numpy':np.__version__,'pandas':pd.__version__},'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'elapsed_seconds':time.monotonic()-START,'production_raw_body_reads':0,'reference_body_reads':0,'production_mutex_operations':0,'real_worker_or_fit_operations':0,'full35_trait_adapter_supplied':False,'execution_admitted':False}
with (R/'core_column_reader_component_controls_receipt_v3_1.json').open('x') as f:json.dump(receipt,f,indent=2);f.write('\n')
print(json.dumps({key:receipt[key] for key in ['status','check_count','case_count','peak_rss_bytes','elapsed_seconds']}))
