#!/usr/bin/env python3
"""Additional small adapter92 identity/timing cases; inherits prior method fixtures."""
import contextlib
import importlib.util
import json
from pathlib import Path
import types

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('controls',HERE/'validation_collector_adapter_checker_v1_1.py')
c=importlib.util.module_from_spec(spec);spec.loader.exec_module(c)
c.FIX=HERE/'validation_collector_adapter_extra_fixtures_v1'
def main():
 c.assert_(not c.FIX.exists(),'fresh_extra_fixture_namespace');c.FIX.mkdir()
 c.sys.path.insert(0,str(c.S));oldurlopen=c.urllib.request.urlopen
 def no_network(*args,**kwargs):raise AssertionError('NETWORK_FORBIDDEN')
 c.urllib.request.urlopen=no_network
 try:
  for p,h in c.EXPECTED.items():c.assert_(c.sha(p)==h,'before_'+p)
  c.reject('prior_output',lambda f:c.write(f['member']['new_munged'],'previous'),expected='PRIOR_PRIVATE_VALIDATION_OUTPUT_PRESERVED_NO_RETRY',promoted=True)
  def sym(f):
   target=Path(f['member']['new_qc']);target.parent.mkdir(parents=True);target.symlink_to(f['body'])
  c.reject('prior_output_symlink',sym,expected='PRIOR_PRIVATE_VALIDATION_OUTPUT_PRESERVED_NO_RETRY')
  def badchecksum(f):
   f['row']['replication_checksum']='md5:'+'0'*32
   with Path(f['plan']['queue']).open(newline='') as h:rows=list(c.csv.DictReader(h,delimiter='\t'))
   rows[0]=f['row'];c.write_tsv(Path(f['plan']['queue']),rows);f['plan']['dependencies_sha256'][f['plan']['queue']]=c.sha(f['plan']['queue'])
  c.reject('inner_MD5_mismatch',badchecksum,expected='MD5 mismatch')
  def missing_fields(rows):return [{k:v for k,v in r.items() if k!='sebeta'} for r in rows]
  f=c.fixture('missing_required_header',transform=missing_fields);err,_=c.adapter_run(f)
  c.assert_(err is not None and 'lacks required fields' in err,'missing_required_header_reject')
  c.assert_(not Path(f['member']['new_munged']).exists(),'missing_required_header_no_promotion')
  c.cases.append(dict(name='missing_required_header',status='PASS_REJECTION',error=err))
  f=c.fixture('zero_survivors',transform=lambda rows:[{**r,'rsids':'rs99999'} for r in rows]);err,_=c.adapter_run(f)
  c.assert_(err is not None and 'too few HM3 variants' in err,'zero_survivors_reject')
  c.assert_(not Path(f['member']['new_munged']).exists(),'zero_survivors_output_removed')
  c.cases.append(dict(name='zero_survivors',status='PASS_REJECTION',error=err))
  f=c.fixture('inner_EOF_before_promotion');events=[]
  def observe(a,f):
   origload=a.load
   def load(name,path):
    m=origload(name,path)
    if name=='_original_validation_streaming':
     orig=m.open_verified_gzip_text
     @contextlib.contextmanager
     def reader(**kw):
      with orig(**kw) as (text,receipt):yield text,receipt
      c.assert_(receipt['verification_status']=='PASS','actual_inner_EOF_PASS')
      c.assert_(receipt['observed_sha256']==f['member']['original_source_identity']['sha256'],'actual_inner_EOF_full_SHA')
      c.assert_(not Path(f['member']['new_munged']).exists(),'actual_inner_EOF_precedes_promotion')
      events.append('verified_EOF_before_promotion')
     m.open_verified_gzip_text=reader
    return m
   a.load=load
  err,_=c.adapter_run(f,observe);c.assert_(err is None and events==['verified_EOF_before_promotion'],'healthy_EOF_timing')
  c.cases.append(dict(name='inner_EOF_before_promotion',status='PASS',events=events,error=err))
  f=c.fixture('receipt_after_parse_before_hash')
  def parse_swap(a,f):
   loads=json.loads
   def load_swap(text,*args,**kw):
    parsed=loads(text,*args,**kw)
    if isinstance(parsed,dict) and parsed.get('pipeline_status')=='STREAM_HARMONIZE_DIRECT_MUNGE_PASS':
     changed={**parsed,'pipeline_status':'SYNTHETIC_CHANGED_AFTER_PARSE'}
     Path(f['member']['new_receipt']).write_text(json.dumps(changed,indent=2)+'\n')
    return parsed
   a.json=types.SimpleNamespace(loads=load_swap,dumps=json.dumps)
  err,_=c.adapter_run(f,parse_swap)
  c.assert_(err is None,'receipt_parse_before_hash_FALSE_ACCEPT')
  final=json.loads(Path(f['member']['new_receipt']).read_text())
  c.assert_(final['pipeline_status']=='SYNTHETIC_CHANGED_AFTER_PARSE','consumed_original_receipt_status_not_bound')
  side=json.loads(Path(f['member']['adapter_receipt']).read_text())
  c.assert_(side['output_sha256'][f['member']['new_receipt']]==c.sha(f['member']['new_receipt']),'sidecar_credits_changed_postparse_receipt')
  c.cases.append(dict(name='receipt_after_parse_before_hash',status='FALSE_ACCEPT',error=err))
  for p,h in c.EXPECTED.items():c.assert_(c.sha(p)==h,'after_'+p)
 finally:c.urllib.request.urlopen=oldurlopen
 result=dict(scope='TINY_ONLY_ADAPTER_BOUNDARY',case_count=len(c.cases),check_count=len(c.checks),cases=c.cases,checks=c.checks,
  source_sha256=c.EXPECTED,runtime=dict(python=c.sys.version,executable=c.sys.executable),
  all_expectations_pass=all(x['passed'] for x in c.checks),no_network=True,no_production_bodies=True,no_mutex=True,
  status='METHOD_EOF_CONTROLS_PASS_CONSUMED_RECEIPT_IDENTITY_FALSE_ACCEPT_MEASURED')
 c.save(HERE/'validation_collector_adapter_extra_controls_v1.json',result)
 print(json.dumps({k:result[k] for k in ['status','case_count','check_count']}))
if __name__=='__main__':main()
