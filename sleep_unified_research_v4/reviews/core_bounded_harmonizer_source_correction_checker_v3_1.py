"""Narrow source-binding correction fixtures; reader4 eligibility separate."""
import ast
import copy
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path
import resource
import sys
import time

import pandas as pd
import numpy as np

P=Path(__file__).resolve().parents[1];R=P/'reviews';S=P/'scripts';sys.path.insert(0,str(S))
OLD=S/'core_bounded_harmonizer_v1.py';NEW=S/'core_bounded_harmonizer_v3.py';READER=S/'core_column_reader_v4.py';ORIGINAL=P.parent/'scripts/01_harmonize.py'
ROOT=R/'core_bounded_harmonizer_source_correction_controls_v3_1';ROOT.mkdir(exist_ok=False)
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest();CHECKS=[];CASES=[];START=time.monotonic()
def check(name,condition,detail=None):
    CHECKS.append({'control':name,'pass':bool(condition),'detail':detail})
    if not condition:raise AssertionError(name)
BINDS=[OLD,NEW,READER,ORIGINAL,R/'core_bounded_harmonizer_critical_review_seal_v1.json']
BEFORE={str(p):sha(p) for p in BINDS}
check('v3_exact_code_identity',BEFORE[str(NEW)]=='3b87e849e3001a2a2d30fc0971dd3bd97847736abb9bed0dff521d1b5e1d0210')
check('reader4_exact_binding_pending_review',BEFORE[str(READER)]=='1e5a7f86ab51cd9cec2d705144ca7623bd8e57e218c3889cf1a94b4db7c5b375')
check('declared_runtime',sys.version_info[:3]==(3,11,11) and pd.__version__=='2.2.3' and np.__version__=='1.26.4')
seal=json.loads((R/'core_bounded_harmonizer_critical_review_seal_v1.json').read_text())
for path,value in seal['artifacts'].items():check('preserved_v1_artifact_'+Path(path).name,sha(path)==value['sha256'] and Path(path).stat().st_size==value['bytes'])
oldtree=ast.parse(OLD.read_text());newtree=ast.parse(NEW.read_text())
oldnodes={x.name:x for x in oldtree.body if isinstance(x,(ast.FunctionDef,ast.ClassDef))};newnodes={x.name:x for x in newtree.body if isinstance(x,(ast.FunctionDef,ast.ClassDef))}
for name,node in oldnodes.items():
    if name!='harmonize':check('unchanged_global_scientific_AST_'+name,ast.dump(node,include_attributes=False)==ast.dump(newnodes[name],include_attributes=False))
harmonize=copy.deepcopy(newnodes['harmonize']);harmonize.body=[x for i,x in enumerate(harmonize.body) if i not in [1,2,3,4,5,15]]
check('harmonize_otherwise_exact_AST_after_source_correction',ast.dump(harmonize,include_attributes=False)==ast.dump(oldnodes['harmonize'],include_attributes=False))
spec=importlib.util.spec_from_file_location('source_corrected_candidate',NEW);c=importlib.util.module_from_spec(spec);spec.loader.exec_module(c)
import core_column_reader_v4 as reader
loader=c.load_original;digest=c.digest;reader_digest=reader.digest;read=pd.read_csv
module,_=loader(ORIGINAL);sniff=module.sniff_separator
config=ROOT/'analysis_panel.tsv';config.write_text('trait_id\ttype\tbuild\tn_total\tsource_note\nx\tcontinuous\thg19\t1000\tOWN_SYNTHETIC_FIXTURE\n')
sealed='SNP\tCHR\tBP\tA1\tA2\tbeta\tb\tSE\tP\tN\nrs1\t1\t42\ta\tc\t1\t9\t1\t0.5\t1000\n'
unsealed='SNP\tCHR\tBP\tA1\tA2\tb\tSE\tP\tN\nrs1\t1\t42\ta\tc\t9\t1\t0.5\t1000\n'
source=ROOT/'source.tsv';source.write_text(sealed);expected=sha(source)
base=['--trait','x','--config',str(config)]
oldargv=sys.argv;sys.argv=[str(ORIGINAL),*base,'--infile',str(source),'--outdir',str(ROOT/'original')]
try:module.main()
finally:sys.argv=oldargv
oldoutput=gzip.decompress((ROOT/'original/x.harmonized.tsv.gz').read_bytes());oldqc=(ROOT/'original/x.qc.txt').read_bytes()

def execute(name,infile_args,wanted=expected,healthy=False,mutate_header=False,error=None):
    folder=ROOT/name;folder.mkdir();events=[]
    def trace_digest(path):
        result=digest(path);events.append({'operation':'digest','path':str(path),'sha256':result});return result
    def trace_read(path,*args,**kwargs):
        result=read(path,*args,**kwargs);events.append({'operation':'read_csv','path':str(path),'header':kwargs.get('nrows')==0})
        if mutate_header and str(path)==str(source) and kwargs.get('nrows')==0:source.write_text(sealed)
        return result
    def trace_load(path):
        events.append({'operation':'load_original','path':str(path)});m,b=loader(path)
        old_sniff=m.sniff_separator
        def trace_sniff(path):events.append({'operation':'sniff','path':str(path)});return old_sniff(path)
        m.sniff_separator=trace_sniff;return m,b
    c.digest=trace_digest;reader.digest=trace_digest;c.load_original=trace_load;pd.read_csv=trace_read
    try:
        try:receipt=c.harmonize(ORIGINAL,[*base,*infile_args,'--outdir',str(folder/'candidate')],folder/'spool',1,wanted);status='SUCCESS';failure=None
        except (Exception,SystemExit) as e:receipt=None;status='FAILED';failure=type(e).__name__+': '+str(e)
    finally:c.digest=digest;reader.digest=reader_digest;c.load_original=loader;pd.read_csv=read
    if healthy:
        check(name+'_healthy_success',status=='SUCCESS',failure)
        newoutput=gzip.decompress((folder/'candidate/x.harmonized.tsv.gz').read_bytes());newqc=(folder/'candidate/x.qc.txt').read_bytes()
        check(name+'_original_output_QC_exact',newoutput==oldoutput and newqc==oldqc)
        firstread=next(i for i,x in enumerate(events) if x['operation'] in ['sniff','read_csv'] and x['path']==str(source))
        check(name+'_correct_digest_precedes_first_raw_source_read',events[0]['operation']=='digest' and events[0]['path']==str(source) and events[0]['sha256']==wanted and firstread>0)
        check(name+'_all_actual_source_read_paths_bound',all(x['path']==str(source) for x in events if x['operation']=='sniff' or x['operation']=='read_csv' and x['path']!=str(config)))
    else:
        check(name+'_fails_closed',status=='FAILED' and error in failure,failure)
        check(name+'_no_candidate_output_or_manifest',not (folder/'candidate').exists() and not (folder/'spool/global_harmonization_candidate_receipt.json').exists())
        if error!='PARSED_ORIGINAL_INFILE_DIFFERS':check(name+'_rejects_before_module_or_source_read',all(x['operation']=='digest' for x in events))
    CASES.append({'name':name,'infile_arguments':infile_args,'expected_sha256':wanted,'status':status,'failure':failure,'events':events,'reader4_eligibility':'PENDING_SEPARATE_INDEPENDENT_COMPONENT_SEAL'})
    return events

execute('healthy_separate_tokens',['--infile',str(source)],healthy=True)
execute('healthy_equals_form',['--infile='+str(source)],healthy=True)
source.write_text(unsealed)
events=execute('former_prefix_header_false_PASS',['--infile',str(source)],mutate_header=True,error='IDENTITY_DIFFERS_BEFORE_SCHEMA_PREFIX')
check('former_witness_now_no_header_or_prefix_read',not any(x['operation'] in ['sniff','read_csv','load_original'] for x in events))
source.write_text(sealed)
execute('wrong_expected_SHA',['--infile',str(source)],wanted='0'*64,error='IDENTITY_DIFFERS_BEFORE_SCHEMA_PREFIX')
for name,args,error in [('missing',[],'EXACT_SINGLE'),('missing_value',['--infile'],'INFILE_VALUE_REQUIRED'),('empty_equals',['--infile='],'EXACT_SINGLE'),('empty_value',['--infile',''],'EXACT_SINGLE'),('duplicate',['--infile',str(source),'--infile',str(source)],'EXACT_SINGLE'),('mixed_duplicate',['--infile='+str(source),'--infile',str(source)],'EXACT_SINGLE'),('abbreviation_only',['--infil',str(source)],'EXACT_SINGLE')]:execute(name,args,error=error)
link=ROOT/'same_bytes_symlink.tsv';link.symlink_to(source)
execute('symlink_final_component',['--infile',str(link)],error='IDENTITY_DIFFERS_BEFORE_SCHEMA_PREFIX')
execute('directory_not_file',['--infile',str(ROOT)],error='IDENTITY_DIFFERS_BEFORE_SCHEMA_PREFIX')
execute('missing_file',['--infile',str(ROOT/'missing.tsv')],error='IDENTITY_DIFFERS_BEFORE_SCHEMA_PREFIX')

# Original argparse accepts abbreviations. An exact token plus abbreviated
# override must reject parsed-source mismatch before column/full-data/output.
other=ROOT/'other_source.tsv';other.write_text(unsealed)
events=execute('parsed_source_abbreviation_override',['--infile',str(source),'--infil',str(other)],error='PARSED_ORIGINAL_INFILE_DIFFERS')
check('parsed_mismatch_no_full_source_parse',not any(x['operation']=='read_csv' and x['path']==str(other) and not x['header'] for x in events))
check('parsed_mismatch_interface_qualification_measured',any(x['operation']=='sniff' and x['path']==str(other) for x in events))
AFTER={p:sha(p) for p in BEFORE};check('all_bound_sources_and_original_rejection_seal_unchanged',BEFORE==AFTER)
receipt={'status':'PASS_V3_SOURCE_CORRECTION_ONLY_READER4_COMPONENT_ELIGIBILITY_PENDING','check_count':len(CHECKS),'checks':CHECKS,'cases':CASES,'sources_before':BEFORE,'sources_after':AFTER,'inherited_v1_global_methods_checks':105,'inherited_v1_global_methods_cases':22,'reader4_source_sha256':BEFORE[str(READER)],'reader4_component_review_status':'PENDING_SEPARATE_INDEPENDENT_REVIEW','interface_qualification':'Exact separate-token/equals routes hash source before sniff/header. Original argparse abbreviated override can read another header before parsed-source gate rejects; no full parse/output/manifest occurs. Frozen production commands must use the exact reviewed route with no abbreviated overrides.','peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'elapsed_seconds':time.monotonic()-START,'production_body_reads':0,'reference_body_reads':0,'production_mutex_or_workers_or_fits':0,'execution_admitted':False}
with (R/'core_bounded_harmonizer_source_correction_controls_receipt_v3_1.json').open('x') as f:json.dump(receipt,f,indent=2);f.write('\n')
print(json.dumps({'status':receipt['status'],'check_count':len(CHECKS),'case_count':len(CASES),'peak_rss_bytes':receipt['peak_rss_bytes']}))
