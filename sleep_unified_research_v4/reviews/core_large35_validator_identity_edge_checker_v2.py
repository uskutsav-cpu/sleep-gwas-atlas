"""Actual validator79 on invented tiny files; no production reads or workers."""
import contextlib
import gzip
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import sys

P=Path(__file__).resolve().parents[1];R=P/'reviews';S=P/'scripts'
sys.path.insert(0,str(S));sys.dont_write_bytecode=True
ROOT=R/'core_large35_validator_identity_edge_controls_v2';ROOT.mkdir(exist_ok=False)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,v):Path(p).write_text(json.dumps(v,indent=2)+'\n')
def loadmod():
 spec=importlib.util.spec_from_file_location('validator_identity_edge79',S/'79_validate_core_large35_bounded_replay_v2.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
CASES=[]
for case in ['healthy','QC_after_parse_before_capture','harmonized_after_comparison_before_capture','munged_after_comparison_before_capture','QC_after_capture']:
 root=ROOT/case;root.mkdir();m=loadmod();raw=root/'raw_identity_metadata';raw.write_text('INVENTED-NOT-GWAS\n');h=root/'new.harmonized.gz';oldh=root/'old.harmonized.gz';z=root/'new.munged.gz';oldz=root/'old.munged.gz'
 hb=b'SNP\tCHR\tBP\tA1\tA2\tFRQ\tBETA\tSE\tP\tN\nOWN_METADATA\t1\t1\ta\tc\t0.2\t0.1\t0.01\t0.1\t100\n';zb=b'SNP\tA1\tA2\tZ\tN\nOWN_METADATA\ta\tc\t1\t100\n'
 for p,b in [(h,hb),(oldh,hb),(z,zb),(oldz,zb)]:
  with gzip.open(p,'wb') as f:f.write(b)
 qc=root/'QC.txt';qc.write_text('rows_in\t1\nrows_out\t1\nstep\tdropped\tremaining\nvalid\t0\t1\n');gate=root/'source_gate.json';out=root/'comparison.json'
 member={'trait_id':'OWN_METADATA','harmonized':str(h),'munged':str(z),'harmonization_qc':str(qc),'source_gate_receipt':str(gate),'original_design':{'raw':{'resolved_path':str(raw),'bytes':raw.stat().st_size,'sealed_verified_sha256':sha(raw)},'prefilter':None,'harmonized':{'path':str(oldh),'sealed_sha256':sha(oldh)},'munged':{'path':str(oldz),'sealed_sha256':sha(oldz)},'harmonized_rows':1,'historical_munge_log':{'output_total_rows':1},'historical_QC':{'metadata':{'rows_in':'1','rows_out':'1'},'ordered_steps':[{'reason':'valid','dropped':0,'remaining':1}]}}}
 plan=root/'plan.json';save(plan,{'members':[member],'dependencies_sha256':{}});ph=sha(plan);save(gate,{'status':'EXACT_CORE_RAW_SOURCE_GATE_PASS','trait':'OWN_METADATA','plan_sha256':ph,'raw_sha256':sha(raw)})
 calls=[0];normal_qc=m.qc;normal_h=m.compare_harmonized;normal_z=m.compare_munged
 def mutated_qc(path):
  observed=normal_qc(path)
  if case=='QC_after_parse_before_capture':qc.write_text('rows_in\t1\nrows_out\t999\nstep\tdropped\tremaining\nvalid\t0\t999\n')
  return observed
 def mutated_h(*a):
  v=normal_h(*a)
  if case=='harmonized_after_comparison_before_capture':
   with gzip.open(h,'wb') as f:f.write(hb.replace(b'\t0.1\t0.01',b'\t9.1\t0.01'))
  return v
 def mutated_z(*a):
  v=normal_z(*a)
  if case=='munged_after_comparison_before_capture':
   with gzip.open(z,'wb') as f:f.write(zb.replace(b'\t1\t100',b'\t9\t100'))
  return v
 def bindings(p):
  calls[0]+=1
  if case=='QC_after_capture' and calls[0]==2:qc.write_text('rows_in\t1\nrows_out\t999\nstep\tdropped\tremaining\nvalid\t0\t999\n')
 m.qc=mutated_qc;m.compare_harmonized=mutated_h;m.compare_munged=mutated_z;m.check_bindings=bindings
 previous=sys.argv;sys.argv=['validator79','--plan',str(plan),'--plan-sha',ph,'--trait','OWN_METADATA','--mode','compare','--out',str(out)];error=None
 try:
  with contextlib.redirect_stdout(io.StringIO()):m.main()
 except BaseException as e:error=type(e).__name__+': '+str(e)
 finally:sys.argv=previous
 result=json.loads(out.read_text());now=normal_qc(qc);match=now==({'rows_in':'1','rows_out':'1'},[{'reason':'valid','dropped':0,'remaining':1}]);accepted=result['status']=='EXACT_CORE_HARMONIZED_MUNGED_CONTENT_AND_QC_REPLAY_PASS'
 expected=case in ['healthy','QC_after_parse_before_capture'];assert accepted==expected,(case,error,result)
 CASES.append({'case':case,'validator_accepts':accepted,'current_QC_matches_original':match,'error':error,'comparison_receipt':str(out),'comparison_receipt_sha256':sha(out),'observed_identity_map':result.get('output_identity_sha256'),'is_false_PASS_witness':accepted and not match})
receipt={'status':'MEASURED_VALIDATOR79_QC_CONSUMPTION_IDENTITY_ESCAPE','source_sha256':sha(S/'79_validate_core_large35_bounded_replay_v2.py'),'cases':CASES,'actual_workers':0,'production_bodies':0,'real_mutex':0,'fits':0,'execution_admitted':False}
with (R/'core_large35_validator_identity_edge_receipt_v2.json').open('x') as f:json.dump(receipt,f,indent=2);f.write('\n')
print(json.dumps(receipt,indent=2))
