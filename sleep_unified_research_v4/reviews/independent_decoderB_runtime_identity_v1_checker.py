#!/usr/bin/env python3
"""Independently verify current runtime identity; never inspect scientific data.

Own scandir physical traversal does not follow directory links. Each runtime
regular byte is independently hashed; link targets must first exactly match
the fixed contained profile and are proved by the full physical file hashes.
No producer census function is imported or used.
"""
from collections import Counter
from datetime import datetime,timezone
import hashlib
import json
import os
from pathlib import Path
import re
import resource
import stat
import subprocess
import sys
import time
sys.dont_write_bytecode=True
P=Path(__file__).resolve().parents[1];R=P/'reviews'
PROFILE=P/'source_provenance/signed_reference_decoderB_runtime_identity_v1.json'
PPIN='e1e9755f90c306aad9294b6066168ceab4e366b30038df99fb24b963d5bfea6c'
RECEIPT=P/'logs/signed_reference_decoderB_runtime_identity_audit_receipt_v2.json'
RPIN='577ed1c4638b0d6f8a592ac8bdf1c9c91d63311e7cdd0a8bbd29c8aaebb15a40'
STOP=P/'logs/signed_reference_decoderB_runtime_identity_audit_failure_v1.json'
SPIN='52d9e297e89917e8829d618df5fecca61d03875667133a521d56bd4bcbad7e34'
OUT=R/'independent_decoderB_runtime_identity_v1_controls.json'
START=time.monotonic();CHECKS=[];PRESERVED={};CENSUSES=[]
def check(v,label):
 CHECKS.append({'label':label,'passed':bool(v)})
 if not v:raise AssertionError(label)
def budget():
 if time.monotonic()-START>300:raise RuntimeError('OWN_REVIEW_300SEC_CAP')
 if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss>128*2**20:raise RuntimeError('OWN_REVIEW_128MiB_RSS_CAP')
def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for block in iter(lambda:f.read(65536),b''):h.update(block)
 return h.hexdigest()
def fixed(path,expected):
 path=Path(path);check(path.is_file() and not path.is_symlink(),'regular_fixed:'+str(path));h=sha(path);check(h==expected,'fixed_sha:'+str(path));PRESERVED[str(path)]=h;return h
def load(path,expected):
 fixed(path,expected);r=json.loads(Path(path).read_text());check(sha(path)==expected,'fixed_postparse:'+str(path));return r
def canonical(value):return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
def inside(name,target):
 # Component membership, not lexical prefix: /a/tree_other is excluded.
 return Path(name).is_relative_to(Path(target))
def directory_tree(target,regular,file_links,dir_links,dirs):
 selected={name:value for name,value in regular.items() if inside(name,target)}
 value={'regular_files':selected,'physical_directories':sorted(name for name in dirs if inside(name,target)),
        'file_symlinks':{name:{k:r[k] for k in ['literal_link','resolved_target','resolved_sha256']} for name,r in file_links.items() if inside(name,target)},
        'directory_symlinks':{name:{k:r[k] for k in ['literal_link','resolved_target']} for name,r in dir_links.items() if inside(name,target)}}
 return canonical(value),len(selected),sum(r['bytes'] for r in selected.values())
def synthetic_controls():
 base='/OWN_INVENTED_RUNTIME';target=base+'/tree';a=target+'/x';outside=base+'/tree_other/x'
 reg={outside:{'bytes':11,'sha256':'b'*64},a:{'bytes':3,'sha256':'a'*64}}
 links={target+'/link':{'literal_link':'x','resolved_target':a,'resolved_sha256':'a'*64,'target_is_file':True}}
 dl={target+'/sub_link':{'literal_link':'sub','resolved_target':target+'/sub','target_is_directory':True,'UNCONSUMED_EXTRA':True}}
 dirs=[base+'/tree_other',target+'/sub',base,target]
 expected={'regular_files':{a:reg[a]},'physical_directories':[target,target+'/sub'],'file_symlinks':{target+'/link':{'literal_link':'x','resolved_target':a,'resolved_sha256':'a'*64}},'directory_symlinks':{target+'/sub_link':{'literal_link':'sub','resolved_target':target+'/sub'}}}
 value,n,b=directory_tree(target,reg,links,dl,dirs);check(value==canonical(expected) and n==1 and b==3,'invented_exact_tree_schema_and_component_boundary')
 check(directory_tree(target,dict(reversed(list(reg.items()))),links,dl,list(reversed(dirs)))==(value,n,b),'invented_sorted_compact_canonical_order')
 reg[outside]['sha256']='f'*64;check(directory_tree(target,reg,links,dl,dirs)==(value,n,b),'invented_outside_descendant_excluded')
 links[target+'/link']['literal_link']='./x';check(directory_tree(target,reg,links,dl,dirs)[0]!=value,'invented_literal_link_changes_tree_sha');links[target+'/link']['literal_link']='x'
 reg[a]['sha256']='c'*64;check(directory_tree(target,reg,links,dl,dirs)[0]!=value,'invented_regular_sha_changes_tree_sha');reg[a]['sha256']='a'*64
 dl[target+'/sub_link']['resolved_target']=target;check(directory_tree(target,reg,links,dl,dirs)[0]!=value,'invented_directory_link_target_changes_tree_sha')
def census(profile,label):
 start=time.monotonic();root=Path(profile['environment_root']);check(root.is_dir() and not root.is_symlink() and not any(p.is_symlink() for p in root.parents),label+':regular_root_and_parents')
 regular=profile['regular_files'];fl=profile['symlinks'];dl=profile['symlink_directories'];seen_r=set();seen_fl=set();seen_dl=set();dirs=set();pending=[root];facts={};total=0;nodes=0
 while pending:
  folder=pending.pop();dirs.add(str(folder))
  with os.scandir(folder) as entries:
   for entry in entries:
    path=Path(entry.path);name=str(path);before=entry.stat(follow_symlinks=False);kind=before.st_mode;nodes+=1
    if stat.S_ISDIR(kind):pending.append(path)
    elif stat.S_ISREG(kind):
     if name not in regular:raise RuntimeError('NEW_RUNTIME_REGULAR_FILE:'+name)
     expected=regular[name]
     if set(expected)!= {'bytes','sha256'} or type(expected['bytes']) is not int or expected['bytes']<0 or not re.fullmatch('[0-9a-f]{64}',expected['sha256']):raise RuntimeError('FIXED_REGULAR_SCHEMA_INVALID:'+name)
     if before.st_size!=expected['bytes'] or sha(path)!=expected['sha256']:raise RuntimeError('CURRENT_REGULAR_RUNTIME_IDENTITY_DIFFERS:'+name)
     after=path.lstat();fact=(after.st_dev,after.st_ino,after.st_size,after.st_mtime_ns,after.st_ctime_ns)
     if fact!=(before.st_dev,before.st_ino,before.st_size,before.st_mtime_ns,before.st_ctime_ns) or not stat.S_ISREG(after.st_mode):raise RuntimeError('REGULAR_RUNTIME_CHANGED_WHILE_HASHING:'+name)
     seen_r.add(name);facts[name]=fact;total+=after.st_size
     if total>2**30:raise RuntimeError('OWN_RUNTIME_ONE_GiB_CENSUS_LIMIT')
    elif stat.S_ISLNK(kind):
     literal=os.readlink(path);target=path.resolve(strict=True)
     # Reject before reading a link target outside the fixed runtime identity.
     if not target.is_relative_to(root):raise RuntimeError('CURRENT_LINK_OUTSIDE_RUNTIME_NOT_READ:'+name)
     if target.is_dir():
      expected=dl.get(name)
      if expected is None or expected['literal_link']!=literal or expected['resolved_target']!=str(target) or expected['target_is_directory'] is not True:raise RuntimeError('CURRENT_DIRECTORY_LITERAL_OR_TARGET_DIFFERS:'+name)
      seen_dl.add(name)
     elif target.is_file() and not target.is_symlink():
      expected=fl.get(name)
      if expected is None or expected['literal_link']!=literal or expected['resolved_target']!=str(target) or expected['target_is_file'] is not True or str(target) not in regular or regular[str(target)]['sha256']!=expected['resolved_sha256']:raise RuntimeError('CURRENT_FILE_LITERAL_OR_TARGET_DIFFERS:'+name)
      seen_fl.add(name)
     else:raise RuntimeError('UNSUPPORTED_RESOLVED_RUNTIME_LINK:'+name)
     after=path.lstat()
     if (before.st_dev,before.st_ino,before.st_mtime_ns,before.st_ctime_ns)!=(after.st_dev,after.st_ino,after.st_mtime_ns,after.st_ctime_ns) or os.readlink(path)!=literal or path.resolve(strict=True)!=target:raise RuntimeError('RUNTIME_LINK_CHANGED_DURING_CENSUS:'+name)
    else:raise RuntimeError('UNSUPPORTED_CURRENT_RUNTIME_NODE:'+name)
    if nodes%256==0:budget()
 check(seen_r==set(regular) and seen_fl==set(fl) and seen_dl==set(dl),label+':complete_exact_file_and_link_membership')
 check(dirs==set(profile['physical_directories'])==set(profile['directories']),label+':complete_exact_physical_directory_membership')
 check(len(seen_r)==26334 and total==390617547 and len(seen_fl)==1225 and len(seen_dl)==2,label+':exact_current_counts_bytes')
 # Link target bytes were independently hashed as physical regular files in
 # this same pass. Recheck metadata to avoid accepting a target swap afterward.
 for value in fl.values():
  name=value['resolved_target'];s=Path(name).lstat()
  if name not in seen_r or (s.st_dev,s.st_ino,s.st_size,s.st_mtime_ns,s.st_ctime_ns)!=facts[name]:raise RuntimeError('LINK_TARGET_CHANGED_AFTER_OWN_BYTE_HASH:'+name)
 tree_results={}
 for name,value in dl.items():
  target=value['resolved_target'];check(target in dirs,label+':directory_link_physical_target:'+name)
  digest,count,size=directory_tree(target,regular,fl,dl,dirs)
  check(digest==value['tree_sha256'] and count==value['physical_regular_file_count'] and size==value['physical_regular_file_bytes'],label+':full_contained_tree_sha_count_bytes:'+name)
  tree_results[name]={'literal_link':value['literal_link'],'resolved_target':target,'tree_sha256':digest,'regular_file_count':count,'regular_file_bytes':size}
 exact={'regular_files':regular,'file_symlinks':fl,'physical_directories':sorted(dirs),'directory_symlinks':dl}
 result={'pass':label,'regular_file_count':len(seen_r),'regular_file_bytes':total,'file_symlink_count':len(seen_fl),'directory_symlink_count':len(seen_dl),'physical_directory_count':len(dirs),'canonical_complete_census_sha256':canonical(exact),'directory_trees':tree_results,'elapsed_seconds':time.monotonic()-start,'runtime_regular_bytes_independently_hashed':total,'file_link_target_sha_verified_through_current_physical_hashes':True}
 print(json.dumps({'progress':'CENSUS_COMPLETE','pass':label,'files':len(seen_r),'bytes':total,'elapsed_seconds':result['elapsed_seconds']}),flush=True);budget();return result
def main():
 check(not OUT.exists(),'fresh_independent_receipt');profile=load(PROFILE,PPIN);receipt=load(RECEIPT,RPIN);stop=load(STOP,SPIN)
 fixed(receipt['producer_script'],receipt['producer_script_sha256']);fixed(stop['producer_script'],stop['producer_script_sha256'])
 check(receipt['profile']==str(PROFILE) and receipt['profile_sha256']==PPIN and receipt['preserved_v1_stop_receipt_sha256']==SPIN,'producer_receipt_binds_exact_profile_and_preserved_stop')
 check(receipt['complete_two_censuses_equal'] is True and receipt['versions_before_after_equal'] is True and receipt['source_stage_execution_admitted'] is False and receipt['independent_runtime_profile_review_pass'] is False,'producer_claim_exact_scope_not_independent_or_source_admission')
 check(stop['profile_or_success_receipt_written'] is False and stop['runtime_modified_installed_or_copied'] is False and stop['status']=='STOPPED_BEFORE_PROFILE_WRITE_EXISTING_DIRECTORY_LINK_UNSUPPORTED','v1_stopped_and_preserved_qualified')
 check(profile['complete_environment_census'] is True and profile['all_before_after_identities_unchanged'] is True and profile['versions']==profile['versions_before']==profile['versions_after'],'profile_complete_two_equal_producer_censuses_and_versions')
 check(profile['regular_file_count']==26334 and profile['regular_file_bytes']==390617547 and profile['symlink_count']==1225 and profile['symlink_directory_count']==2,'profile_declared_cardinality')
 check(profile['historical_per_trait_runtime_attestation_recovered'] is False and profile['exact_historical_binary_recovery_claimed'] is False and profile['data_body_reads']==0 and profile['worker_fit_calls']==0 and profile['runtime_installs_or_copies']==0,'historical_and_science_claims_bounded')
 synthetic_controls();CENSUSES.append(census(profile,'independent_before_query'))
 root=Path(profile['environment_root']);logical=root/'bin/python';binary=root/'bin/python3.9'
 check(logical.is_symlink() and os.readlink(logical)=='python3.9' and logical.resolve(strict=True)==binary and not binary.is_symlink(),'exact_logical_python_link_and_resolved_regular_binary')
 check(sha(binary)==profile['inherited_resolved_python_binary_sha256']=='7dffb088cd3027e48f0127ced6f206e06111abed3b9457e580b6ebbf593c10ba','current_python_binary_sha')
 # Separately written version query; -B and explicit single-thread variables.
 query="import json,sys,numpy,pandas,bitarray;print(json.dumps({'python':sys.version,'executable':sys.executable,'prefix':sys.prefix,'base_prefix':sys.base_prefix,'numpy':numpy.__version__,'pandas':pandas.__version__,'bitarray':bitarray.__version__}))"
 command=[str(logical),'-B','-c',query];environment={**os.environ,'PYTHONDONTWRITEBYTECODE':'1','OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','MKL_NUM_THREADS':'1','NUMEXPR_NUM_THREADS':'1','VECLIB_MAXIMUM_THREADS':'1'}
 run=subprocess.run(command,env=environment,capture_output=True,text=True,check=True,timeout=60);versions=json.loads(run.stdout)
 check(versions==profile['versions'],'independently_queried_exact_versions_and_logical_environment')
 check(versions['python'].startswith('3.9.23') and (versions['numpy'],versions['pandas'],versions['bitarray'])==('1.21.5','1.3.3','2.8.3'),'declared3_9_23_1_21_5_1_3_3_2_8_3')
 check(versions['executable']==str(logical) and versions['prefix']==str(root),'exact_logical_interpreter_prefix')
 CENSUSES.append(census(profile,'independent_after_query'));check({k:v for k,v in CENSUSES[0].items() if k not in ['pass','elapsed_seconds']}=={k:v for k,v in CENSUSES[1].items() if k not in ['pass','elapsed_seconds']},'two_independent_current_censuses_equal_and_match_profile')
 for p,h in PRESERVED.items():check(sha(p)==h,'final_preserved_identity:'+p)
 budget();record={'schema':'independent_complete_current_decoderB_runtime_identity_review','verdict':'QUALIFIED_CURRENT_RUNTIME_IDENTITY_PASS_ONLY','check_count':len(CHECKS),'checks':CHECKS,'censuses':CENSUSES,'two_independent_full_current_censuses_equal':True,'profile_sha256':PPIN,'producer_success_receipt_sha256':RPIN,'preserved_v1_stop_sha256':SPIN,'preserved_producer_profile_receipt_sha256':PRESERVED,'version_query_command':command,'version_query_exact_result':versions,'version_query_stderr':run.stderr,'version_query_BLAS_threads':1,'version_query_pycache_disabled':True,'native_python_binary_sha256':sha(binary),'elapsed_seconds':time.monotonic()-START,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'checker_sha256':sha(__file__),'source_genotype_reference_GWAS_body_reads':0,'scientific_decoder_calls':0,'estimator_calls':0,'runtime_installs_copies_modifications':0,'actual_mutex':0,'network':0,'source_execution_admission':False,'historical_runtime_execution_attestation':False,'exact_historical_runtime_recovery':False,'system_library_identity_coverage':False,'scientific_decoder_qualification':False,'source_eligibility_qualification':False,'qualification':'Current physical environment, contained literal symlinks/trees, declared versions and binary identities only. System/dynamic libraries outside environment, historical executions, scientific decoding/orientation and source rights/eligibility remain independent gates. Current equality cannot certify an earlier runtime epoch. Link targets were proved through freshly hashed current physical regular identities without following any external target.'}
 with OUT.open('x') as f:json.dump(record,f,indent=2,allow_nan=False);f.write('\n');f.flush();os.fsync(f.fileno())
 print(json.dumps({'receipt_sha256':sha(OUT),'assertions':len(CHECKS),'elapsed_seconds':record['elapsed_seconds'],'peak_rss_bytes':record['peak_rss_bytes']}),flush=True)
if __name__=='__main__':main()
