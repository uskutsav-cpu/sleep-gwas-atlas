"""Independent bounded fixtures; never accepts production raw inputs."""
import ast
from contextlib import contextmanager
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
import warnings

import numpy as np
import pandas as pd

P=Path(__file__).resolve().parents[1]
REPO=P.parent
SOURCE=P/'scripts/core_column_reader_v1.py'
ORIGINAL=REPO/'scripts/01_harmonize.py'
ROOT=P/'reviews/core_column_reader_component_controls_v1'
ROOT.mkdir(exist_ok=False)
CHECKS=[]
CASES=[]
START=time.monotonic()

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

BEFORE={str(x):sha(x) for x in [SOURCE,ORIGINAL,P/'reviews/core_bounded_parser_methods_review_v1.md',P/'reviews/core_bounded_parser_methods_review_v1.json',P/'reviews/core_bounded_parser_methods_review_seal_v1.json']}
assert BEFORE[str(SOURCE)]=='0b9ac982d1537ccee7f8572b98286d0ab00cc769ecbd401627f41cf6c0961e5e'
spec=importlib.util.spec_from_file_location('candidate_core_column_reader',SOURCE)
candidate=importlib.util.module_from_spec(spec)
spec.loader.exec_module(candidate)
assert sys.version_info[:3]==(3,11,11) and pd.__version__=='2.2.3' and np.__version__=='1.26.4'

# This oracle executes unchanged original whole-DataFrame normalization AST,
# after the original full selected-column parse and canonical rename.
tree=ast.parse(ORIGINAL.read_text())
main=next(x for x in tree.body if isinstance(x,ast.FunctionDef) and x.name=='main')
nodes=[x for x in main.body if 430<=x.lineno<=454]
ORACLE=compile(ast.fix_missing_locations(ast.Module(body=nodes,type_ignores=[])),str(ORIGINAL),'exec')

def check(name,condition,detail=None):
    CHECKS.append({'control':name,'pass':bool(condition),'detail':detail})
    if not condition:raise AssertionError(name)

def fingerprint(s):
    h=hashlib.sha256()
    h.update(str(s.dtype).encode())
    h.update(np.asarray(s.isna(),dtype=np.uint8).tobytes())
    if isinstance(s.dtype,np.dtype) and s.dtype.kind in 'biufc':
        h.update(s.to_numpy(copy=False).tobytes())
    else:
        for x in s:
            h.update(type(x).__module__.encode()+b'.'+type(x).__name__.encode()+b':')
            if isinstance(x,(float,np.floating)):
                h.update(struct.pack('>d',float(x)))
            else:h.update(repr(x).encode())
            h.update(b'\0')
    return h.hexdigest()

def kwargs(engine,sep,columns,skip=0):
    out={'sep':sep,'comment':None,'skiprows':skip,'usecols':list(dict.fromkeys(columns.values()))}
    if engine=='python':out['engine']='python'
    else:out['low_memory']=True
    return out

def original(path,kw,columns):
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')
        raw=pd.read_csv(path,**kw)
        parsed={c:{'dtype':str(raw[c].dtype),'fingerprint':fingerprint(raw[c])} for c in raw}
        builds=[str(x) for x in raw[columns['BUILD']].dropna().unique() if str(x).strip()] if 'BUILD' in columns else []
        out=raw.rename(columns={v:k for k,v in columns.items()})
        env={'out':out,'pd':pd,'np':np,'log':lambda x:None}
        exec(ORACLE,env)
    return env['out'],parsed,builds,[str(x.message) for x in w]

def restore(manifest,target):
    parts=manifest['typed_pieces'][target]
    loaded=[]
    for x in parts:
        check(target+'_piece_hash_'+Path(x['path']).name,sha(x['path'])==x['sha256'])
        with open(x['path'],'rb') as f:s=pickle.load(f)
        check(target+'_piece_dtype_'+Path(x['path']).name,str(s.dtype)==x['dtype'])
        check(target+'_piece_index_'+Path(x['path']).name,s.index.equals(pd.RangeIndex(x['start'],x['start']+x['rows'])))
        loaded.append(s)
    if loaded:return pd.concat(loaded)
    return pd.Series([],dtype=manifest['global_dtype'][target])

def run_case(name,text,columns,engine='c',sep='\t',skip=0,compressed=False):
    case=ROOT/name;case.mkdir()
    path=case/('source.txt.gz' if compressed else 'source.txt')
    content=text.encode()
    if compressed:
        with gzip.open(path,'wb') as f:f.write(content)
    else:path.write_bytes(content)
    kw=kwargs(engine,sep,columns,skip)
    out,parsed,builds,oldwarnings=original(path,kw,columns)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')
        manifest=candidate.freeze_columns(path,kw,columns,case/'spool',len(out))
    check(name+'_source_before_after',manifest['source_sha256_before']==manifest['source_sha256_after']==sha(path))
    check(name+'_BUILD_global_values',manifest['original_BUILD_values']==builds)
    targets={('BETA' if c=='OR' else c) for c in columns}
    record={'name':name,'rows':len(out),'physical_source_bytes':path.stat().st_size,'source_sha256':sha(path),'kwargs':kw,'columns':columns,'parsed_columns':parsed,'original_warnings':oldwarnings,'candidate_warnings':[str(x.message) for x in w],'normalized':{}}
    for target in sorted(targets):
        restored=restore(manifest,target)
        pd.testing.assert_series_equal(out[target],restored,check_names=False,check_exact=True)
        check(name+'_'+target+'_global_dtype',str(out[target].dtype)==manifest['global_dtype'][target])
        check(name+'_'+target+'_exact_fingerprint',fingerprint(out[target])==fingerprint(restored))
        check(name+'_'+target+'_NA_mask',np.array_equal(out[target].isna().to_numpy(),restored.isna().to_numpy()))
        record['normalized'][target]={'dtype':str(out[target].dtype),'fingerprint':fingerprint(out[target]),'NA_count':int(out[target].isna().sum())}
    if engine=='c':
        widths=[];blocks=[]
        for selected in [kw['usecols'],[kw['usecols'][0]]]:
            r=pd.read_csv(path,**dict(kw,usecols=selected,iterator=True))
            widths.append(int(r._engine._reader.table_width))
            a=r._engine._reader.read_low_memory(None)
            blocks.append([len(next(iter(z.values()))) for z in a]);r.close()
        check(name+'_actual_full_table_width',widths[0]==widths[1]==manifest['physical_shape_gate']['full_header_width'])
        check(name+'_actual_native_block_boundaries',blocks[0]==blocks[1])
        record['native_full_width']=widths;record['native_block_rows']=blocks
    CASES.append(record)
    return manifest

BASE={'SNP':'snp','A1':'a1','A2':'a2','CHR':'chr','BP':'bp','BETA':'beta','SE':'se','P':'p','FRQ':'frq','N':'n','BUILD':'build'}
run_case('C_numeric_NA_bounds','snp\ta1\ta2\tchr\tbp\tbeta\tse\tp\tfrq\tn\tbuild\n001\ta\t C \t1\t00042\t-0.0\t1e-2\t0.5\tNA\t18446744073709551615\thg19\n1e0\tT\tG\tchr2\t9223372036854775808\tbad\tNaN\tinf\t0.25\t-1\tGRCh37\nNULL\tNA\tN\t3\t1\t1e-310\t0\t-0\tNaN\tNaN\tNA\n',BASE)
ORCOL=dict(BASE);ORCOL.pop('BETA');ORCOL['OR']='or'
run_case('C_OR_full_log','snp\ta1\ta2\tchr\tbp\tse\tp\tfrq\tn\tbuild\tor\nrs1\ta\tc\t1\t2\t1\t.5\t.2\t10\thg19\t1.0000000000000002\nrs2\tt\tg\t2\t3\t1\t.5\t.2\t10\thg19\t0\nrs3\ta\tc\t3\t4\t1\t.5\t.2\t10\thg38\t-1\nrs4\ta\tc\t4\t5\t1\t.5\t.2\t10\tNA\tNaN\n',ORCOL,compressed=True)
MIN={'SNP':'snp','A1':'a1','A2':'a2','BETA':'beta'}
run_case('C_comma_quotes_metadata','metadata\nmetadata2\nsnp,a1,a2,beta,unused\n" RS,1 "," a ",c,1,"quoted,extra"\nrs2,t,g,-0.0,x\n',MIN,sep=',',skip=2)

# Wide physical input gives native C boundaries much smaller than50k.  Four
# native blocks mix booleans, numeric values, float lexical forms and objects.
width=64;nrows=32783
head=['snp','a1','a2','beta','n','lex']+['unused'+str(i) for i in range(width-6)]
wide=ROOT/'wide_native_source.txt'
with wide.open('w') as f:
    f.write('\t'.join(head)+'\n')
    for i in range(nrows):
        beta='True' if i<8192 else ('1' if i<16384 else ('1.00000000000000001' if i<24576 else 'bad'))
        n='18446744073709551615' if i<16384 else '-1'
        lex='001' if i<24576 else 'late'
        f.write('\t'.join(['rs'+str(i),'a','c',beta,n,lex]+['0']*(width-6))+'\n')
widecols=dict(MIN,N='n',BUILD='lex')
run_case('C_wide_native_inference',wide.read_text(),widecols)

run_case('Python_global_fallback','snp a1 a2 beta n build\n001 a c 001 001 hg19\n1e0 t g 1e0 1e0 NA\nlate a c bad True hg38\nNULL N NA NaN False NULL\n',dict(MIN,N='n',BUILD='build'),engine='python',sep=r'\s+')
run_case('Python_short_blank_quoted','metadata\nsnp a1 a2 beta n\nrs1 a c 1 10\n\n  \nrs2 t g\nrs3 a c 2\n"RS 4" a c 3 4\nrs5 a c 1 2 extra\n',dict(MIN,N='n'),engine='python',sep=r'\s+',skip=1,compressed=True)
run_case('Python_all_NA','snp a1 a2 beta n\nNA NA NA NA NA\nNULL null NaN nan N/A\n',dict(MIN,N='n'),engine='python',sep=r'\s+')
run_case('Python_bool_numeric','snp a1 a2 beta n\nrs1 a c True True\nrs2 t g False False\nrs3 a c 1 0\nrs4 a c bad bad\n',dict(MIN,N='n'),engine='python',sep=r'\s+')
run_case('Python_cross_50k_fallback','snp a1 a2 beta\n'+'001 a c 001\n'*50001+'1e0 t g bad\n',MIN,engine='python',sep=r'\s+')
run_case('Python_empty_get_lines_chunk','snp a1 a2 beta\n'+'\n'*50002+'001 a c 001\n1e0 t g bad\n',MIN,engine='python',sep=r'\s+')
run_case('C_header_only','snp\ta1\ta2\tbeta\n',MIN)
run_case('Python_header_only','snp a1 a2 beta\n',MIN,engine='python',sep=r'\s+')

def expect_reject(name,text,columns=MIN,engine='c',sep='\t',expected=1,transform=None,error_contains=None):
    folder=ROOT/name;folder.mkdir();path=folder/'source.txt'
    path.write_text(text);kw=kwargs(engine,sep,columns)
    if transform:transform(path,kw,folder)
    try:candidate.freeze_columns(path,kw,columns,folder/'spool',expected)
    except Exception as e:
        check(name+'_reject',error_contains is None or error_contains in str(e),type(e).__name__+': '+str(e))
        return
    raise AssertionError(name+'_unexpected_acceptance')

expect_reject('C_ragged','snp\ta1\ta2\tbeta\nrs1\ta\tc\n',error_contains='RAGGED')
expect_reject('C_implicit_index','snp\ta1\ta2\tbeta\n0\trs1\ta\tc\t1\n',error_contains='IMPLICIT_INDEX')
expect_reject('C_multiline','snp\ta1\ta2\tbeta\n"rs\n1"\ta\tc\t1\n',error_contains='MULTILINE')
expect_reject('C_duplicate_header','snp\ta1\ta2\tbeta\tbeta\nrs1\ta\tc\t1\t1\n',error_contains='UNIQUE')
expect_reject('C_expected_rows','snp\ta1\ta2\tbeta\nrs1\ta\tc\t1\n',expected=2,error_contains='COUNT')
expect_reject('Python_expected_rows','snp a1 a2 beta\nrs1 a c 1\n',engine='python',sep=r'\s+',expected=2,error_contains='COUNT')
expect_reject('Python_implicit_index','snp a1 a2 beta\n0 rs1 a c 1\n',engine='python',sep=r'\s+',error_contains='IMPLICIT')
expect_reject('C_BETA_preference','snp\ta1\ta2\tbeta\tor\nrs1\ta\tc\t1\t2\n',columns=dict(MIN,OR='or'),error_contains='BETA_PREFERENCE')
for engine in ['c','python']:
    folder=ROOT/(engine+'_CRC_failure');folder.mkdir();path=folder/'source.txt.gz'
    text='snp\ta1\ta2\tbeta\nrs1\ta\tc\t1\n' if engine=='c' else 'snp a1 a2 beta\nrs1 a c 1\n'
    raw=bytearray(gzip.compress(text.encode()));raw[-8]^=1;path.write_bytes(raw)
    try:candidate.freeze_columns(path,kwargs(engine,'\t' if engine=='c' else r'\s+',MIN),MIN,folder/'spool',1)
    except (gzip.BadGzipFile,EOFError) as e:check(engine+'_CRC_reject',True,str(e))
    else:raise AssertionError(engine+'_CRC_accepted')

# Controlled real file mutation during the shape scan interval exposes whether
# the candidate's recorded initial hash includes its first source operation.
folder=ROOT/'C_shape_scan_source_mutation';folder.mkdir();path=folder/'source.txt'
initial='snp\ta1\ta2\tbeta\nrs1\ta\tc\t1\n';changed=initial.replace('\t1\n','\t9\n')
path.write_text(initial);initial_sha=sha(path)
gate=candidate.c_shape_gate
def mutate_after_shape(*args,**kw):
    result=gate(*args,**kw);path.write_text(changed);return result
candidate.c_shape_gate=mutate_after_shape
try:
    manifest=candidate.freeze_columns(path,kwargs('c','\t',MIN),MIN,folder/'spool',1)
finally:candidate.c_shape_gate=gate
check('C_shape_mutation_false_PASS_witness',manifest['source_sha256_before']==manifest['source_sha256_after']==sha(path) and initial_sha!=sha(path))
mutation={'status':'CANDIDATE_FALSE_PASS_INITIAL_SHAPE_SCAN_MUTATION','initial_sha256':initial_sha,'accepted_manifest_source_before':manifest['source_sha256_before'],'accepted_source_after':manifest['source_sha256_after'],'accepted_BETA':float(restore(manifest,'BETA').iloc[0]),'minimal_fix':'Compute initial digest before any source shape/parser scan, retain final digest, and bind caller expected raw identity.'}
(folder/'witness.json').write_text(json.dumps(mutation,indent=2)+'\n')

folder=ROOT/'C_between_column_mutation';folder.mkdir();path=folder/'source.txt';path.write_text(initial)
read=candidate.pd.read_csv;calls=0
def mutate_between_reads(*args,**kw):
    global calls
    if args and Path(args[0])==path:
        calls+=1
        if calls==2:path.write_text(changed)
    return read(*args,**kw)
candidate.pd.read_csv=mutate_between_reads
try:
    try:candidate.freeze_columns(path,kwargs('c','\t',MIN),MIN,folder/'spool',1)
    except RuntimeError as e:check('C_between_column_mutation_reject','RAW_SOURCE_CHANGED' in str(e),str(e))
    else:raise AssertionError('C_between_column_mutation_accepted')
finally:candidate.pd.read_csv=read

AFTER={str(x):sha(x) for x in map(Path,BEFORE)}
check('all_bound_source_metadata_unchanged',BEFORE==AFTER)
receipt={'status':'QUALIFIED_VALUE_PARITY_WITH_ONE_SOURCE_INTEGRITY_BLOCKER','versions':{'python':sys.version,'numpy':np.__version__,'pandas':pd.__version__,'executable':sys.executable},'source_metadata_before':BEFORE,'source_metadata_after':AFTER,'checks':CHECKS,'check_count':len(CHECKS),'cases':CASES,'false_PASS_source_integrity_witness':mutation,'candidate_rows_per_typed_piece':candidate.ROWS,'fixture_only':True,'production_raw_body_reads':0,'reference_body_reads':0,'real_worker_or_fit_operations':0,'production_mutex_operations':0,'component_only_execution_admitted':False,'elapsed_seconds':time.monotonic()-START,'peak_reviewer_RSS_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
(P/'reviews/core_column_reader_component_controls_receipt_v1.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'status':receipt['status'],'check_count':len(CHECKS),'case_count':len(CASES),'peak_RSS_bytes':receipt['peak_reviewer_RSS_bytes'],'elapsed_seconds':receipt['elapsed_seconds']}))
