#!/usr/bin/env python3
"""Independent tiny-only original47 adapter boundary fixtures. Never production IO."""
import ast
import contextlib
import csv
import gzip
import hashlib
import importlib.util
import io
import json
import os
from pathlib import Path
import shutil
import sys
import urllib.request

REPO=Path(__file__).resolve().parents[2]
P=REPO/'sleep_unified_research_v4'
S=P/'scripts'
R=P/'reviews'
FIX=R/'validation_collector_adapter_fixtures_v1_1'
ADAPTER=S/'92_replay_original_validation_collector_v1.py'
ORIGINAL=REPO/'discovery_extension/scripts/47_stream_replication_sources.py'
STREAM=REPO/'discovery_extension/scripts/streaming_io.py'
TERMINAL=S/'terminal_commit_common_v2.py'
EXPECTED={str(ADAPTER):'15253ace8479334198f2bff06f33a9a5afd690f0c47d3f0ff31a4ad9f372e4fc',
 str(ORIGINAL):'43d00fdf535a2e323f3cd1b9e93f49a29a1683452e67a40fed9e8916e2ffe4a5',
 str(STREAM):'28af3379171194ace67c2c09ffe8f3a46e83a54da1875be72a69535fbb566ac8',
 str(TERMINAL):'9ef14eda37d07b6a3442181820f8909bedd7764aae01fbcaafc966791dd5a7dd'}
checks=[];cases=[]
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def assert_(condition,label):
 checks.append(dict(label=label,passed=bool(condition)))
 if not condition:raise AssertionError(label)
def load(name,path):
 spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
def write(path,content):
 path=Path(path);path.parent.mkdir(parents=True,exist_ok=True);path.write_text(content);return path
def save(path,value):return write(path,json.dumps(value,indent=2,allow_nan=False)+'\n')
def write_tsv(path,rows):
 path.parent.mkdir(parents=True,exist_ok=True)
 with path.open('w',newline='') as f:
  writer=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(rows)
def gzip_table(path,rows):
 text=io.StringIO();writer=csv.DictWriter(text,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(rows)
 path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(gzip.compress(text.getvalue().encode(),compresslevel=1,mtime=0));return path
def adapted_source():
 text=ORIGINAL.read_text();before=ast.parse(text);after=ast.parse(text)
 replacements=[]
 for node in ast.walk(after):
  if isinstance(node,ast.Constant) and type(node.value)==int and node.value in (900_000,700_000):
   replacements.append(node.value);node.value=1
 assert_(sorted(replacements)==[700_000,900_000],'exactly_two_fixture_only_threshold_constants')
 # Revert our two nodes and prove the complete AST unchanged.
 reverted=ast.parse(text.replace('900_000','1').replace('700_000','1'))
 assert_(ast.dump(after,include_attributes=False)==ast.dump(reverted,include_attributes=False),'only_threshold_AST_delta')
 return text.replace('900_000','1').replace('700_000','1')

def finn_rows():
 def row(rs,**kw):
  r=dict(ref='C',alt='A',rsids=rs,beta='0.3',sebeta='0.1',af_alt='0.2');r.update(kw);return r
 rows=[row('rs1'),row('rs1',beta='100'),row('rs2',ref='A',alt='C',beta='0.4',sebeta='0.2'),
  row('rs3',ref='G',alt='T',beta='0.6',sebeta='0.2'),row('rs4',ref='T',alt='G',beta='0.8',sebeta='0.2'),
  row('rs5'),row('rs6'),row('rs7'),row('rs8'),row('rs999'),row('rs40'),
  row('rs9',af_alt='0.01'),row('rs10',af_alt='0.99'),row('rs11',af_alt='-0.1'),row('rs12',af_alt='nan'),
  row('rs13',beta='nan'),row('rs14',sebeta='0'),row('rs15',sebeta='-1'),row('rs16',sebeta='nan'),
  row('rs17',alt='A',ref='T'),row('rs18',alt='N'),row('rs19',alt='A',ref='A'),
  row('rs20',beta='1e308',sebeta='1e-308'),row('rs21',af_alt='bad'),
  row('rs22',beta='bad'),row('rs23',sebeta='bad'),row('rs24',ref='G',alt='A'),
  row('rs25',af_alt='0.010000000001'),row('rs26',af_alt='0.989999999999'),
  row('rs999;rs27;rs28',ref='G',alt='A'), # First matched rs27 rejected; no rs28 alternative.
  row('rs29',sebeta='0'),row('rs29',beta='0.9',sebeta='0.3'), # first invalid duplicate does not claim ID.
  row('prefixrs30suffix',alt='a',ref='c'),row('rs31',af_alt='1.1')]
 return rows
def mvp_rows():
 def row(rs,**kw):
  r=dict(effect_allele='A',other_allele='C',odds_ratio='2',effect_allele_frequency='0.2',p_value='0.5',rsid=rs,r2='0.95',ci_upper='4',ci_lower='1');r.update(kw);return r
 rows=[row('rs1'),row('rs1'),row('rs2',effect_allele='C',other_allele='A'),
  row('rs3',effect_allele='T',other_allele='G'),row('rs4',effect_allele='G',other_allele='T'),
  row('rs5'),row('rs6'),row('rs7'),row('rs8'),row('rs999'),row('rs40'),
  row('rs9',r2='0.9'),row('rs10',r2='0.900000001'),row('rs11',r2='nan'),
  row('rs12',odds_ratio='0'),row('rs13',odds_ratio='nan'),row('rs14',p_value='-1'),
  row('rs15',p_value='1.1'),row('rs16',p_value='nan'),row('rs17',p_value='0'),row('rs18',p_value='1'),
  row('rs19',ci_lower='0'),row('rs20',ci_upper='1'),row('rs21',ci_upper='0.5'),row('rs22',ci_lower='nan'),
  row('rs23',effect_allele_frequency='0.01'),row('rs24',effect_allele_frequency='0.99'),
  row('rs25',effect_allele_frequency='0.010000000001'),row('rs26',effect_allele_frequency='nan'),
  row('rs27',effect_allele='A',other_allele='T'),row('rs28',effect_allele='N'),
  row('rs29',effect_allele='C',other_allele='C'),row('rs30',ci_upper='1.0000000000000002',odds_ratio='1e308'),
  row('rs31',odds_ratio='bad'),row('rs32',r2='bad'),row('rs33',ci_upper='bad'),
  row('rs34',effect_allele_frequency='bad'),row('rs35',effect_allele='A',other_allele='G'),
  row('rs36',odds_ratio='-1'),row('rs36',odds_ratio='1'),row('prefixrs37suffix',effect_allele='a',other_allele='c')]
 return rows

def fixture(name,kind='finn',transform=None):
 root=FIX/name;root.mkdir(parents=True);workspace=root/'workspace';workspace.mkdir()
 science=root/'science';collector=write(science/'47_stream_replication_sources.py',adapted_source());stream=write(science/'streaming_io.py',STREAM.read_text())
 ref=root/'reference.tsv.gz';alleles=root/'alleles.tsv';records=[];ars=[]
 for i in range(1,51):
  chr_,bp=(6,25_000_000) if i==5 else (6,34_000_000) if i==6 else (6,24_999_999) if i==7 else (6,34_000_001) if i==8 else (1,i)
  records.append(dict(SNP='rs'+str(i),CHR=chr_,BP=bp));ars.append(dict(SNP='rs'+str(i),A1='A',A2='T' if i==40 else 'C'))
 gzip_table(ref,records);write_tsv(alleles,ars)
 sid='finngen_r13_fixture_0' if kind=='finn' else 'gwas_catalog_fixture_0'
 rows=finn_rows() if kind=='finn' else mvp_rows()
 if transform:rows=transform(rows)
 body=gzip_table(root/'body.tsv.gz',rows)
 def identity():
  b=body.read_bytes();return dict(bytes=len(b),md5=hashlib.md5(b).hexdigest(),sha256=hashlib.sha256(b).hexdigest())
 orig=identity();row=dict(replication_source_id=sid,source_curation_status='COMPLETE_BEFORE_RESULTS',
  replication_source_url='https://fixture.invalid/'+sid,replication_content_length_bytes=str(orig['bytes']),
  replication_etag=orig['md5'],replication_storage_mode='VERSIONED_REMOTE_STREAMING' if kind=='finn' else 'PUBLIC_STUDY_ACCESSION_STREAMING',
  replication_source_generation='123456',replication_checksum='md5:'+orig['md5'],cases='123',controls='987',
  replication_munged_path='discovery_extension/data/replication/munged/'+sid+'.sumstats.gz',
  replication_receipt_path='discovery_extension/provenance/replication_streaming_receipts/'+sid+'.json')
 qrows=[row]+[{**row,'replication_source_id':'finngen_r13_fixture_'+str(i),'replication_source_url':'https://fixture.invalid/unused_'+str(i)} for i in range(1,13)]
 queue=root/'queue.tsv';write_tsv(queue,qrows)
 # Real terminal helper on only synthetic immutable metadata: no transfer/mutex.
 acq=root/'acquisition';acq.mkdir();individual=[]
 for i in range(26):individual.append(save(acq/(str(i)+'.json'),dict(synthetic_only=True,index=i)))
 primary=save(acq/'primary.json',dict(status='SYNTHETIC_COMPLETED_NO_BODY_ACQUISITION'))
 pending=acq/'PENDING';seal=acq/'seal.json';binding=dict(synthetic_only=True,fixture=name)
 terminal=load('fixture_terminal_'+name,TERMINAL);tc=terminal.TerminalCommit(pending,seal,binding)
 assert_(tc.commit({str(primary):sha(primary)},identity_gate=lambda:None,resource_gate=lambda:None,termination_gate=lambda:None),'synthetic_terminal_commit_'+name)
 member=dict(source_id=sid,body_path=str(body),original_source_identity=orig,frozen_source_metadata=row,
  new_munged=str(workspace/row['replication_munged_path']),new_receipt=str(workspace/row['replication_receipt_path']),
  new_qc=str(workspace/('discovery_extension/results/replication/qc/'+sid+'.tsv')),adapter_receipt=str(workspace/'adapter.json'))
 member['collector_argv']=[str(collector),'--source-id',sid,'--queue',str(queue),'--reference',str(ref),'--hm3-alleles',str(alleles),'--acknowledge-network-gib',format(orig['bytes']/(1<<30),'.17g')]
 deps={str(p):sha(p) for p in [collector,stream,ref,alleles,queue,ADAPTER,TERMINAL]}
 contract=dict(metadata_sha256={str(p):sha(p) for p in [*individual,primary,seal]},individual_receipt_paths=list(map(str,individual)),
  pending_path=str(pending),seal_path=str(seal),binding=binding,primary_receipt_sha256={str(primary):sha(primary)})
 plan=dict(member_count=13,new_estimator_fits=0,members=[member]+[dict(source_id='finngen_r13_fixture_'+str(i)) for i in range(1,13)],
  dependencies_sha256=deps,acquisition_terminal=contract,workspace=str(workspace),original_collector=str(collector),
  original_streaming_io=str(stream),reference=str(ref),alleles=str(alleles),queue=str(queue),acquisition_plan_sha256='f'*64)
 pp=save(root/'plan.json',plan)
 return dict(root=root,workspace=workspace,plan=plan,member=member,row=row,body=body,pp=pp,identity=identity)

def adapter_run(f,hook=None):
 a=load('adapter_'+f['root'].name,ADAPTER)
 if hook:hook(a,f)
 with contextlib.redirect_stdout(io.StringIO()) as out:
  try:a.run(f['pp'],sha(f['pp']),f['member']['source_id']);error=None
  except BaseException as e:error=type(e).__name__+': '+str(e)
 return error,out.getvalue()
def update(f):save(f['pp'],f['plan'])
def qc(path):
 with Path(path).open(newline='') as h:return list(csv.DictReader(h,delimiter='\t'))
def oracle(f):
 root=f['root']/'oracle';root.mkdir();row=f['row'];pp=f['plan'];cp=root/'discovery_extension/scripts/47_stream_replication_sources.py';sp=cp.parent/'streaming_io.py'
 write(cp,Path(pp['original_collector']).read_text());write(sp,STREAM.read_text())
 stream=load('oracle_stream_'+f['root'].name,sp);prior=sys.modules.get('streaming_io');sys.modules['streaming_io']=stream
 try:m=load('oracle_'+f['root'].name,cp)
 finally:
  if prior is None:sys.modules.pop('streaming_io',None)
  else:sys.modules['streaming_io']=prior
 m.ROOT=root
 m.head=lambda u:{'etag':row['replication_etag'],'content-length':row['replication_content_length_bytes'],'x-goog-generation':row['replication_source_generation']}
 @contextlib.contextmanager
 def local(**kwargs):
  with stream.open_verified_gzip_text(local_path=f['body'],expected_md5=kwargs['expected_md5'],expected_size_bytes=kwargs['expected_size_bytes'],timeout_seconds=kwargs['timeout_seconds']) as x:yield x
 m.open_verified_gzip_text=local;argv=f['member']['collector_argv'].copy();argv[0]=str(cp)
 oldargv=sys.argv;cwd=Path.cwd();sys.argv=argv;os.chdir(root)
 try:
  with contextlib.redirect_stdout(io.StringIO()):m.main()
 finally:sys.argv=oldargv;os.chdir(cwd)
 return root

def healthy(name,kind):
 f=fixture(name,kind);err,log=adapter_run(f);assert_(err is None,name+'_adapter_pass')
 baseline=oracle(f);mp=f['member'];out=Path(mp['new_munged']);base=baseline/f['row']['replication_munged_path']
 literal=gzip.decompress(out.read_bytes());assert_(literal==gzip.decompress(base.read_bytes()),name+'_exact_original_ordered_literals')
 oq=qc(baseline/('discovery_extension/results/replication/qc/'+mp['source_id']+'.tsv'));aq=qc(mp['new_qc'])
 assert_([r for r in aq if r['metric']!='munged_output_sha256']==[r for r in oq if r['metric']!='munged_output_sha256'],name+'_exact_original_QC_except_transport_SHA')
 assert_(next(r['value'] for r in aq if r['metric']=='munged_output_sha256')==sha(out),name+'_actual_compressed_QC_SHA')
 expected_n=format(4.0/(1.0/123+1.0/987),'.12g');parsed=list(csv.DictReader(io.StringIO(literal.decode()),delimiter='\t'))
 assert_(all(r['N']==expected_n for r in parsed),name+'_independent_harmonic_N')
 assert_([r['SNP'] for r in parsed][:6]==['rs1','rs2','rs3','rs4','rs7','rs8'],name+'_independent_order_MHC_boundaries')
 if kind=='finn':
  assert_([r['Z'] for r in parsed][:4]==['3','-2','3','-4'],name+'_independent_all_orientation_Z')
  assert_('rs27' not in [r['SNP'] for r in parsed] and 'rs28' not in [r['SNP'] for r in parsed],name+'_first_matched_rsid')
  assert_(next(r['Z'] for r in parsed if r['SNP']=='rs29')=='3',name+'_first_surviving_duplicate')
  assert_('rs25' in [r['SNP'] for r in parsed] and 'rs26' in [r['SNP'] for r in parsed],name+'_strict_MAF_neighbor')
 else:
  assert_([r['Z'] for r in parsed][:4]==['1.95996398454','-1.95996398454','1.95996398454','-1.95996398454'],name+'_independent_OR_CI_Z')
  assert_('rs9' not in [r['SNP'] for r in parsed] and 'rs10' in [r['SNP'] for r in parsed],name+'_strict_R2_neighbor')
  assert_('rs17' in [r['SNP'] for r in parsed] and 'rs18' in [r['SNP'] for r in parsed],name+'_inclusive_P_boundaries')
  assert_(next(r['Z'] for r in parsed if r['SNP']=='rs36')=='0',name+'_OR_one_zero_Z_and_first_survivor')
 side=json.loads(Path(mp['adapter_receipt']).read_text());rec=json.loads(Path(mp['new_receipt']).read_text())
 assert_(rec['full_resolution_local_retention'] is False and side['full_local_source_retained'] is True,name+'_legacy_retention_truth_qualification')
 assert_(rec['source_verification']['http_etag']=='NA' and side['HTTP_HEAD_not_requeried_in_collector'] is True,name+'_no_fake_current_HTTP')
 assert_(set(rec['pipeline_code_sha256'])=={f['plan']['original_streaming_io'],f['plan']['original_collector']},name+'_routed_code_keys_recorded_absolute')
 assert_(all(sha(path)==h for path,h in side['output_sha256'].items()),name+'_sidecar_output_binding')
 cases.append(dict(name=name,kind='original_main_method_comparison',status='PASS',rows=len(parsed),error=err))

def reject(name,mutate=None,hook=None,expected=None,promoted=False):
 f=fixture(name)
 if mutate:mutate(f);update(f)
 err,log=adapter_run(f,hook);assert_(err is not None,name+'_reject')
 if expected:assert_(expected in err,name+'_exact_rejection')
 if not promoted:assert_(not Path(f['member']['new_munged']).exists(),name+'_no_promoted_output')
 assert_(not Path(f['member']['adapter_receipt']).exists(),name+'_no_success_sidecar')
 cases.append(dict(name=name,kind='negative_boundary',status='PASS_REJECTION',error=err))

def main():
 assert_(not FIX.exists(),'fresh_fixture_namespace');FIX.mkdir()
 for p,h in EXPECTED.items():assert_(sha(p)==h,'pinned_before_'+p)
 sys.path.insert(0,str(S));oldurlopen=urllib.request.urlopen
 def no_network(*args,**kwargs):raise AssertionError('NETWORK_FORBIDDEN_IN_TINY_FIXTURES')
 urllib.request.urlopen=no_network
 try:
  healthy('healthy_finn','finn');healthy('healthy_mvp','mvp')
  reject('wrong_original_SHA',lambda f:f['member']['original_source_identity'].__setitem__('sha256','0'*64),expected='FULL_ORIGINAL_VALIDATION_BODY_CHANGED')
  reject('wrong_original_MD5',lambda f:f['member']['original_source_identity'].__setitem__('md5','0'*32),expected='FULL_ORIGINAL_VALIDATION_BODY_CHANGED')
  reject('wrong_original_size',lambda f:f['member']['original_source_identity'].__setitem__('bytes',f['member']['original_source_identity']['bytes']+1),expected='FULL_ORIGINAL_VALIDATION_BODY_CHANGED')
  reject('wrong_dependency',lambda f:f['plan']['dependencies_sha256'].__setitem__(f['plan']['queue'],'0'*64),expected='FROZEN_VALIDATION_INPUT_CHANGED')
  reject('wrong_command',lambda f:f['member']['collector_argv'].append('--all'),expected='FROZEN_COLLECTOR_ARGUMENTS_CHANGED')
  reject('wrong_route',lambda f:f['member'].__setitem__('new_qc',str(f['root']/'outside.tsv')),expected='EXACT_PRIVATE_ORIGINAL_QC_ROUTE_REQUIRED')
  reject('pending_producer',lambda f:write(f['plan']['acquisition_terminal']['pending_path'],'PENDING'),expected='STAGE_TERMINAL_COMMIT_NOT_COMPLETE')
  reject('failed_individual',lambda f:save(str(f['plan']['acquisition_terminal']['individual_receipt_paths'][0])+'.failure.json',dict(failed=True)),expected='SOURCE_FAILURE_VETOES_REPLAY')
  reject('wrong_receipt_count',lambda f:f['plan']['acquisition_terminal'].__setitem__('individual_receipt_paths',f['plan']['acquisition_terminal']['individual_receipt_paths'][:-1]),expected='ALL_THIRTEEN_DOUBLE_RECEIPT_PATHS_REQUIRED')
  def corrupt(f,truncate=False):
   b=bytearray(f['body'].read_bytes())
   if truncate:b=b[:-5]
   else:b[-8]^=1
   f['body'].write_bytes(b);ident=f['identity']();f['member']['original_source_identity']=ident;f['row']['replication_content_length_bytes']=str(ident['bytes']);f['row']['replication_checksum']='md5:'+ident['md5'];f['row']['replication_etag']=ident['md5']
   with Path(f['plan']['queue']).open(newline='') as h:rs=list(csv.DictReader(h,delimiter='\t'))
   rs[0]=f['row'];write_tsv(Path(f['plan']['queue']),rs);f['plan']['dependencies_sha256'][f['plan']['queue']]=sha(f['plan']['queue']);f['member']['collector_argv'][-1]=format(ident['bytes']/(1<<30),'.17g')
  reject('bad_CRC',corrupt,expected='CRC check failed')
  reject('truncated_gzip',lambda f:corrupt(f,True),expected='end-of-stream marker')
  def patch_receipt(a,f):
   original_load=a.load
   def patched_load(name,path):
    module=original_load(name,path)
    if name=='_original_validation_streaming':
     original_open=module.open_verified_gzip_text
     @contextlib.contextmanager
     def wrong_sha(**kwargs):
      with original_open(**kwargs) as (text,receipt):yield text,receipt
      receipt['observed_sha256']='0'*64
     module.open_verified_gzip_text=wrong_sha
    return module
   a.load=patched_load
  reject('inner_EOF_SHA_mismatch',hook=patch_receipt,expected='FULL_LOCAL_GZIP_EOF_IDENTITY_CHANGED')
  def mutate_after_load(a,f):
   original_load=a.load
   def patched(name,path):
    m=original_load(name,path)
    if name=='_original_validation_streaming':f['body'].write_bytes(f['body'].read_bytes()+b'extra')
    return m
   a.load=patched
  reject('source_after_initial_hash',hook=mutate_after_load)
  # Deliberate counterexamples below are measurements, not admitted input plans.
  f=fixture('unregistered_sidecar_route');outside=f['root']/'outside_workspace_adapter.json';f['member']['adapter_receipt']=str(outside);update(f)
  err,_=adapter_run(f);assert_(err is None and outside.exists(),'unregistered_sidecar_ROUTE_FALSE_ACCEPT_WITNESS')
  cases.append(dict(name='unregistered_sidecar_route',kind='measured_gap',status='FALSE_ACCEPT',error=err,outside_private_workspace_file=str(outside)))
  for target in ['new_qc','plan']:
   f=fixture('post_sidecar_'+target);mp=f['member'];priorPath=Path
   def hook(a,f):
    targetpath=Path(mp['new_qc']) if target=='new_qc' else f['pp']
    class SidecarWriter:
     def __init__(self,handle):self.handle=handle
     def __enter__(self):return self
     def write(self,payload):
      value=self.handle.write(payload)
      if target=='new_qc':targetpath.write_text(targetpath.read_text().replace('input_rows\t34','input_rows\t999'))
      else:targetpath.write_text(targetpath.read_text()+' ')
      return value
     def __exit__(self,*args):return self.handle.__exit__(*args)
    class ScopedPath(type(Path())):
     def open(self,*args,**kwargs):
      h=super().open(*args,**kwargs)
      if str(self)==mp['adapter_receipt'] and args and args[0]=='x':return SidecarWriter(h)
      return h
    a.Path=ScopedPath
   err,_=adapter_run(f,hook);assert_(err is None,'post_sidecar_'+target+'_FALSE_ACCEPT_WITNESS')
   if target=='new_qc':
    side=json.loads(Path(mp['adapter_receipt']).read_text());assert_(side['output_sha256'][mp['new_qc']]!=sha(mp['new_qc']),'post_sidecar_QC_stale_success_identity')
   cases.append(dict(name='post_sidecar_'+target,kind='measured_gap',status='FALSE_ACCEPT',error=err))
  # The double-receipt count alone is not an exact full-family map proof.
  f=fixture('duplicate_26_individual_paths');c=f['plan']['acquisition_terminal'];c['individual_receipt_paths']=[c['individual_receipt_paths'][0]]*26;update(f)
  err,_=adapter_run(f);assert_(err is None,'duplicate_26_receipt_paths_FALSE_ACCEPT_WITNESS')
  cases.append(dict(name='duplicate_26_individual_paths',kind='plan_schema_qualification',status='FALSE_ACCEPT_UNIQUE_MAP_NOT_ENFORCED',error=err))
 finally:urllib.request.urlopen=oldurlopen
 for p,h in EXPECTED.items():assert_(sha(p)==h,'pinned_after_'+p)
 result=dict(scope='TINY_SYNTHETIC_ADAPTER_BOUNDARY_ONLY',runtime=dict(python=sys.version,executable=sys.executable),
  source_sha256=EXPECTED,cases=cases,checks=checks,check_count=len(checks),case_count=len(cases),all_expectations_pass=all(c['passed'] for c in checks),
  fixture_threshold_amendment='Only two in-memory/private-copy integer constants 900000 and700000 become1; production originals unchanged. All method blocks/other literals unchanged.',
  no_network=True,no_production_bodies=True,no_production_workers=True,no_mutex=True,no_new_estimator_fits=True,
  status='REJECT_OPERATIONAL_BOUNDARY_WITH_MEASURED_FALSE_ACCEPTS_METHOD_BRANCHES_MATCH')
 save(R/'validation_collector_adapter_controls_v1_1.json',result)
 print(json.dumps(dict(status=result['status'],checks=len(checks),cases=len(cases))))

if __name__=='__main__':main()
