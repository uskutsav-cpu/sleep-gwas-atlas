"""Additional tiny EOF and mapping guards; no production source arguments."""
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

P=Path(__file__).resolve().parents[1]
SOURCE=P/'scripts/core_column_reader_v1.py'
ROOT=P/'reviews/core_column_reader_component_edge_controls_v1'
ROOT.mkdir(exist_ok=False)
spec=importlib.util.spec_from_file_location('component_edge_candidate',SOURCE)
c=importlib.util.module_from_spec(spec);spec.loader.exec_module(c)
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
before=sha(SOURCE)
assert before=='0b9ac982d1537ccee7f8572b98286d0ab00cc769ecbd401627f41cf6c0961e5e'
checks=[]
cols={'SNP':'snp','A1':'a1','A2':'a2','BETA':'beta'}
for engine in ['c','python']:
    path=ROOT/(engine+'.txt.gz')
    text='snp\ta1\ta2\tbeta\nrs1\ta\tc\t1\n' if engine=='c' else 'snp a1 a2 beta\nrs1 a c 1\n'
    path.write_bytes(gzip.compress(text.encode())[:-2])
    kw={'sep':'\t' if engine=='c' else r'\s+','comment':None,'skiprows':0,'usecols':list(cols.values())}
    kw.update({'engine':'python'} if engine=='python' else {'low_memory':True})
    try:c.freeze_columns(path,kw,cols,ROOT/(engine+'_spool'),1)
    except EOFError as e:checks.append({'control':engine+'_truncated_gzip_EOF','pass':True,'error':str(e)})
    else:raise AssertionError('Truncated gzip accepted')
path=ROOT/'normal.txt';path.write_text('snp\ta1\ta2\tbeta\nrs1\ta\tc\t1\n')
kw={'sep':'\t','comment':None,'skiprows':0,'low_memory':True,'usecols':list(cols.values())[::-1]}
try:c.freeze_columns(path,kw,cols,ROOT/'wrong_order_spool',1)
except RuntimeError as e:
    assert 'SELECTED_COLUMN_ORDER' in str(e)
    checks.append({'control':'selected_column_order_guard','pass':True,'error':str(e)})
else:raise AssertionError('Wrong selected order accepted')
alias=dict(cols,OR='beta');kw['usecols']=list(dict.fromkeys(alias.values()))
try:c.freeze_columns(path,kw,alias,ROOT/'alias_spool',1)
except RuntimeError as e:
    assert 'ONE_TO_ONE' in str(e)
    checks.append({'control':'one_to_one_map_guard','pass':True,'error':str(e)})
else:raise AssertionError('Non-one-to-one mapping accepted')
assert before==sha(SOURCE)
result={'status':'PASS_BOUNDED_EDGE_GUARDS','checks':checks,'check_count':len(checks),'source_sha256_before':before,'source_sha256_after':sha(SOURCE),'production_raw_body_reads':0,'reference_body_reads':0,'real_worker_or_fit_operations':0,'production_mutex_operations':0,'fixture_only':True,'execution_admitted':False,'executable':sys.executable}
(P/'reviews/core_column_reader_component_edge_receipt_v1.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'status':result['status'],'check_count':len(checks)}))
