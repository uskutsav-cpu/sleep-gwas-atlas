#!/usr/bin/env python3
"""Narrow additive intended-seal/validator-persistence controls; invented data only."""
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
from unittest.mock import patch

P=Path(__file__).resolve().parents[1]
path=P/'reviews/independent_original13_pipeline_checker_v2_2.py'
spec=importlib.util.spec_from_file_location('_independent13_component_controls',path)
B=importlib.util.module_from_spec(spec);spec.loader.exec_module(B)
B.ROOT=P.parents[1]/'independent_original13_pipeline_terminal_fixture_v2_1'
out=P/'reviews/independent_original13_pipeline_terminal_controls_v2_1.json'
if B.ROOT.exists() or out.exists():raise RuntimeError('OWN_NEW_NAMESPACE_REQUIRED')
B.ROOT.mkdir()
terminal=sys.modules['terminal_commit_common_v2'];original=terminal.save_new
def altered(path,value):
 original(path,value)
 if Path(path).name=='validation_pipeline_terminal_v2.json':
  # Mutation before Terminal2 captures its own persisted seal digest. The
  # controller must nevertheless preserve its PREWRITE intended identity.
  r=json.loads(Path(path).read_text());r['binding']['executor_sha256']='0'*64
  B.rewrite(path,r)
with patch.object(terminal,'save_new',altered):
 committed,error=B.controller_case('changed_seal_after_write')
 B.check('PREWRITE_intended_seal_rejects_changed_persisted_seal',not committed and error is not None)
for label,target in [('changed_validator_receipt','receipt'),('changed_raw_after_validator_write','body')]:
 d,ssd,namespace,plan,a,m=B.validator_fixture(label);original_sha=B.V.sha;fired=[False]
 def injected(path):
  if str(path)==m['comparison_receipt'] and not fired[0]:
   fired[0]=True;q=Path(path if target=='receipt' else m['body_path']);q.write_bytes(q.read_bytes()+b'INTENTIONAL_OWN_FAULT')
  return original_sha(path)
 with patch.object(B.V,'sha',injected):
  B.expect_fail(label,lambda:B.V.run(a.plan,a.plan_sha,m['source_id'],'compare',Path(m['comparison_receipt'])))
d,ssd,namespace,plan,a=B.fixture('declared_runtime_link')
logical=d/'declared-python-link';logical.symlink_to(Path(plan['python']).resolve())
plan['python']=str(logical);profile=plan['collector_runtime_profile'];profile['logical_executable']=str(logical)
profile['declared_symlinks']=[{'path':str(logical),'target':str(logical.readlink())}]
B.A.runtime_gate(plan);B.check('explicit_known_runtime_link_is_allowed',True)
profile['declared_symlinks'][0]['target']='UNREGISTERED_TARGET'
B.expect_fail('changed_runtime_link_literal_is_rejected',lambda:B.A.runtime_gate(plan))
d=B.ROOT/'scientific_symlink';d.mkdir();regular=d/'regular';regular.write_text('OWN');link=d/'link';link.symlink_to(regular)
B.expect_fail('scientific_symlink_is_rejected',lambda:B.V.regular(link))
rss,size=B.bounded()
receipt={'schema':'independent_original13_pipeline_terminal_controls_v2_1','status':'PASS_SIX_FOCUSED_FAULT_OR_ROUTE_CASES',
 'checks':B.CHECKS,'check_count':len(B.CHECKS),'controller_cases':B.MAP,'synthetic_fixture_bytes':size,'peak_RSS_bytes':rss,
 'candidate_sha256':{str(B.S/k):B.sha(B.S/k) for k in B.EXPECTED},'inherited_checker_sha256':B.sha(path),
 'production_body_reads':0,'network_calls':0,'native_workers':0,'real_heavy_mutex_operations':0,
 'qualification':'Metadata/controller mocks and tiny invented gzip only; no actual prepared processing plan or source producer admission.'}
B.save(out,receipt)
print(json.dumps({'receipt':str(out),'sha256':B.sha(out),'checks':len(B.CHECKS),'peak_RSS_bytes':rss}))
