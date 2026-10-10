"""Compare frozen full original01 to candidate using only owned fixtures."""
import ast
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path
import resource
import sys
import time

import numpy as np
import pandas as pd

P=Path(__file__).resolve().parents[1];R=P/'reviews';S=P/'scripts';REPO=P.parent
sys.path.insert(0,str(S));ORIGINAL=REPO/'scripts/01_harmonize.py';SOURCE=S/'core_bounded_harmonizer_v1.py'
ROOT=R/'core_bounded_harmonizer_methods_controls_v1';ROOT.mkdir(exist_ok=False)
START=time.monotonic();CHECKS=[];CASES=[]
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
column_receipt=json.loads((R/'core_column_reader_component_controls_receipt_v3_1.json').read_text())
BINDS=[SOURCE,ORIGINAL,S/'core_column_reader_v3.py',REPO/'scripts/variant_map.py',REPO/'scripts/liftover_chain.py',R/'core_column_reader_component_review_seal_v3.json']+[Path(p) for p in column_receipt['sources_before']]
BEFORE={str(p):sha(p) for p in BINDS}
assert BEFORE[str(SOURCE)]=='9aa49e2239c52d16cd206b15058000f7e0b5842d6103f920abd1013ebc1703d4'
assert sys.version_info[:3]==(3,11,11) and pd.__version__=='2.2.3' and np.__version__=='1.26.4'
spec=importlib.util.spec_from_file_location('bounded_methods_candidate',SOURCE);c=importlib.util.module_from_spec(spec);spec.loader.exec_module(c)
original,body=c.load_original(ORIGINAL)
def check(name,condition,detail=None):CHECKS.append({'control':name,'pass':bool(condition),'detail':detail})
def table(header,rows,sep='\t'):return sep.join(header)+'\n'+''.join(sep.join(map(str,row))+'\n' for row in rows)
BASE=['SNP','CHR','BP','A1','A2','BETA','SE','P','N']
def row(i=1,n='1000',beta='0.25',p='.5',a1='a',a2='c',chrom=1,bp=None):return ['rs'+str(i),chrom,bp or i,a1,a2,beta,1,p,n]
def invoke(module,argv):
    old=sys.argv;sys.argv=[str(ORIGINAL),*argv]
    try:module.main();return {'status':'SUCCESS'}
    except (Exception,SystemExit) as e:return {'status':'FAILED','error_type':type(e).__name__,'error':str(e)}
    finally:sys.argv=old
def run(name,text,metadata=None,extras=None,schema=None,rows=None,compressed=False,expect_failure=False):
    folder=ROOT/name;folder.mkdir();source=folder/('source.tsv.gz' if compressed else 'source.tsv')
    if compressed:
        with gzip.open(source,'wb') as f:f.write(text.encode())
    else:source.write_text(text)
    config=folder/'analysis_panel.tsv';meta={'trait_id':'x','type':'continuous','build':'hg19','n_total':'1000','ncase':'','ncontrol':'','source_note':'OWN_SYNTHETIC_FIXTURE'};meta.update(metadata or {})
    config.write_text(table(list(meta),[list(meta.values())]))
    if schema:(folder/'gwas_schemas.tsv').write_text(table(['trait_id','sample_size'],[['x',schema]]))
    extra=list(extras or []);argv=['--trait','x','--config',str(config),'--infile',str(source),*extra]
    original_result=invoke(original,[*argv,'--outdir',str(folder/'original')])
    if rows is None:
        sep,_,skip=original.sniff_separator(str(source));kwargs={'sep':sep,'skiprows':skip,'comment':None}
        kwargs.update({'engine':'python'} if sep==r'\s+' else {'low_memory':True});rows=len(pd.read_csv(source,**kwargs))
    try:
        receipt=c.harmonize(ORIGINAL,[*argv,'--outdir',str(folder/'candidate')],folder/'spool',rows,sha(source));candidate_result={'status':'SUCCESS'}
    except (Exception,SystemExit) as e:receipt=None;candidate_result={'status':'FAILED','error_type':type(e).__name__,'error':str(e)}
    check(name+'_same_success_or_rejection',original_result['status']==candidate_result['status'],{'original':original_result,'candidate':candidate_result})
    check(name+'_expected_status',original_result['status']==('FAILED' if expect_failure else 'SUCCESS'))
    record={'name':name,'rows':rows,'source_sha256':sha(source),'original':original_result,'candidate':candidate_result,'original_argv':argv,'receipt':receipt}
    if original_result['status']==candidate_result['status']=='SUCCESS':
        oldout=gzip.decompress((folder/'original/x.harmonized.tsv.gz').read_bytes());newout=gzip.decompress((folder/'candidate/x.harmonized.tsv.gz').read_bytes())
        oldqc=(folder/'original/x.qc.txt').read_bytes();newqc=(folder/'candidate/x.qc.txt').read_bytes()
        check(name+'_exact_decompressed_ordered_output',oldout==newout,{'original_sha256':hashlib.sha256(oldout).hexdigest(),'candidate_sha256':hashlib.sha256(newout).hexdigest()})
        check(name+'_exact_full_ordered_QC_report',oldqc==newqc,{'original_sha256':hashlib.sha256(oldqc).hexdigest(),'candidate_sha256':hashlib.sha256(newqc).hexdigest()})
        check(name+'_all_saved_frame_and_column_manifest_hashes',all(sha(x['path'])==x['sha256'] for x in receipt['post_liftover_frames']+receipt['post_ordinary_QC_frames']) and sha(folder/'spool/columns/column_preparation_manifest.json')==receipt['column_manifest_sha256'])
        record.update(decompressed_output_sha256=hashlib.sha256(newout).hexdigest(),QC_sha256=hashlib.sha256(newqc).hexdigest(),rows_out=receipt['rows_out'])
        if rows<=30:record['decompressed_output']=newout.decode();record['QC_report']=newqc.decode()
    CASES.append(record)
    return record

# All scientific ordinary thresholds are evaluated in the frozen original code.
head=BASE+['INFO','FRQ','OR','BUILD'];good=row(1)+[.9,.25,-1,'hg19']
ordinary=[good,row(2,p=0)+[.9,.25,2,'GRCh37'],row(3,a1='a',a2='t')+[.9,.25,2,'hg19'],row(4,beta='inf')+[.9,.25,2,'hg19'],row(5)+[.6,.25,2,'hg19'],row(6)+[.9,.01,2,'hg19'],row(7,chrom=6,bp=25000000)+[.9,.25,2,'hg19'],row(1,n=1000,beta=9)+[.9,.25,2,'hg19'],row(8,n=499)+[.9,.25,2,'hg19'],row(9,n=500)+[.9,.25,2,'hg19'],row(10,a1='N')+[.9,.25,2,'hg19'],row(11,chrom=23)+[.9,.25,2,'hg19'],row(12)+[1,.99,2,'hg19']]
run('ordinary_thresholds_BETA_preference',table(head,ordinary),compressed=True)
run('Python_whole_global_inference',table(BASE,[row(1,n='001',beta='1.00000000000000001'),row(2,n='1e0',beta='-0.0'),row(3,n='bad',beta='bad'),row(4,n=1000,beta='1e-310')],sep=' '))
run('continuous_all_N_missing_config_fallback',table(BASE,[row(1,n='NA'),row(2,n='NULL')]),metadata={'n_total':'2345.5'})
run('continuous_no_N_config',table(BASE[:-1],[row(1)[:-1],row(2)[:-1]]),metadata={'n_total':'1234'})
run('binary_Neff_precedence',table(BASE+['NEFF','NEFFDIV2','NCASE','NCONTROL'],[row(1)+[1000,1,1,1],row(2)+[499,1000,1000,1000],row(3)+[500,1000,1000,1000]]),metadata={'type':'binary','ncase':'9000','ncontrol':'9000'})
run('binary_half_schema',table(BASE+['Effective_N'],[row(1)+[250.5],row(2)+[500.5],row(3)+[0]]),metadata={'type':'binary'},schema='DERIVED_N_EFF=2*Effective_N')
run('binary_Ncase_Ncontrol_npwhere',table(BASE+['NCASE','NCONTROL'],[row(1)+[10,20],row(2)+[0,20],row(3)+[1000,2000],row(4)+['NA',30]]),metadata={'type':'binary'})
run('binary_config_constant',table(BASE,[row(1),row(2)]),metadata={'type':'binary','ncase':'1000','ncontrol':'2000'})
run('all_ordinary_QC_removed_header_only',table(BASE[:-1],[row(1,beta='bad')[:-1]]),metadata={'n_total':'1000'})
run('global_N_infinite_max_failure',table(BASE,[row(1,n='inf'),row(2,n=1000)]),expect_failure=True)
run('global_BUILD_mixed_failure',table(BASE+['BUILD'],[row(1)+['hg19'],row(2)+['hg38']]),expect_failure=True)
run('all_N_nonpositive_failure',table(BASE,[row(1,n=0),row(2,n=-1)]),expect_failure=True)

for kind in ['duplicate_first_low_N','first_ordinary_invalid','late_global_max','late_global_any','exhausted_first_chunk_global_float']:
    values=[]
    for i in range(50003):
        r=row(i+1,n='NA' if kind=='late_global_any' and i<50002 else 1000)
        if kind=='duplicate_first_low_N':
            if i==0:r=row(1,n=1)
            if i==50002:r=row(1,n=1000,beta=9)
        elif kind=='first_ordinary_invalid':
            if i==0:r=row(1,beta='bad')
            if i==50002:r=row(1,beta=9)
        elif kind=='late_global_max' and i==50002:r[-1]=2001
        elif kind=='exhausted_first_chunk_global_float':
            if i<50000:r[5]='bad'
            if i==50002:r[-1]='1000.5'
        values.append(r)
    record=run('boundary_'+kind,table(BASE,values),metadata={'n_total':'99999'},rows=50003)
    if record.get('receipt'):
        steps=record['receipt']['original_QC_steps']
        if kind=='duplicate_first_low_N':check('duplicate_selected_before_N_then_rejected',record['rows_out']==50001 and any(reason=='duplicate SNP ID' and dropped==1 for reason,dropped,_ in steps))
        if kind=='late_global_max':check('late_max_globally_drops_prior_chunks',record['rows_out']==1)
        if kind=='late_global_any':check('late_any_suppresses_config_fallback_globally',record['rows_out']==1)
    del values,record

# Tiny complete synthetic maps carry exactly the original validated HM3 schema.
def make_map(name,rows):
    path=ROOT/(name+'.tsv.gz')
    with gzip.open(path,'wt') as f:f.write(table(['SNP','CHR','BP','A1','A2'],rows))
    provenance={'schema_version':'atlas.hm3-grch37-variant-map.v1','map_scope':'HAPMAP3_ONLY','genome_build':'GRCh37/hg19','map_sha256':sha(path),'map_bytes':path.stat().st_size,'mapped_rows':len(rows),'fixture_role':'OWN_SYNTHETIC_REFERENCE_ONLY'}
    Path(str(path)+'.provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    return path
reference=make_map('own_map',[['rs1',1,1,'A','C'],['rs2',1,2,'A','C'],['rs3',1,3,'A','C'],['rs4',1,3,'A','C'],['rs5',2,90,'G','T'],['rs6',1,101,'A','C']])
def map_args(strategy,path=reference):return ['--variant-map',str(path),'--variant-map-strategy',strategy,'--expected-variant-map-sha256',sha(path),'--expected-variant-map-bytes',str(path.stat().st_size)]
coord_header=['SNP','A1','A2','BETA','SE','P','N']
run('coordinate_label_case_and_mapping',table(coord_header,[['chr1:1:A:C','a','c',1,1,.5,1000],['CHR1:2:A:C','a','c',1,1,.5,1000],['1_2','a','c',1,1,.5,1000],['1:3','a','c',1,1,.5,1000],['NA','a','c',1,1,.5,1000]]),extras=map_args('BY_COORD_ALLELES'))
run('rsid_map_coordinate_free',table(coord_header,[['RS1','a','c',1,1,.5,1000],['rs2','t','g',1,1,.5,1000],['rs999','a','c',1,1,.5,1000]]),metadata={'build':'unresolved'},extras=map_args('BY_RSID_ALLELES'))
badmap=make_map('required_before_QC_duplicate_map',[['rs1',1,1,'A','C'],['rs7',1,7,'A','C'],['rs7',1,8,'A','C']])
run('map_required_keys_before_ordinary_QC_failure',table(coord_header,[['rs1','a','c',1,1,.5,1000],['rs7','a','c',1,1,0,1000]]),extras=map_args('BY_RSID_ALLELES',badmap),expect_failure=True)
chain=ROOT/'own_chain.gz'
with gzip.open(chain,'wt') as f:f.write('chain 1 chr1 1000 + 0 10 chr1 1000 + 100 110 1\n10\n\nchain 1 chr2 1000 + 0 10 chr2 100 - 10 20 2\n10\n\nchain 1 chr3 1000 + 0 10 chrX 1000 + 100 110 3\n10\n\nchain 1 chr1 1000 + 2 4 chr1 1000 + 200 202 4\n2\n')
life_args=['--liftover-chain',str(chain),'--source-build','hg38','--expected-liftover-chain-sha256',sha(chain),'--expected-liftover-chain-bytes',str(chain.stat().st_size)]
lifecases=[row(1,chrom=1,bp=1),row(2,chrom=2,bp=1,a1='a',a2='c'),row(3,chrom=2,bp=2,a1='a',a2='t'),row(4,chrom=1,bp=3),row(5,chrom=3,bp=1),row(6,chrom=23,bp=1),row(7,chrom=4,bp=1),row(8,chrom=2,bp=3,p=0)]
run('liftover_full_reverse_early_drop_QC',table(BASE+['BUILD'],[x+['hg38'] for x in lifecases]),metadata={'build':'hg38'},extras=life_args)
combined=run('liftover_then_global_map_keys',table(BASE+['BUILD'],[row(6,chrom=1,bp=1)+['hg38'],row(5,chrom=2,bp=1)+['hg38'],row(999,chrom=4,bp=1)+['hg38']]),metadata={'build':'hg38'},extras=life_args+map_args('BY_RSID_ALLELES'))
if combined.get('receipt'):check('required_key_union_after_liftover_drops',combined['receipt']['global_required_map_key_count']==2)

# DiskSet must preserve the original exact set membership and union semantics.
disk=c.DiskSet(ROOT/'own_exact_set.sqlite3');wanted={'rs1','rs2',(1,42,'A','C'),(1,42,'G','T')}
disk.update(wanted);disk.update(wanted)
check('DiskSet_exact_union_cardinality',len(disk)==len(wanted));check('DiskSet_exact_membership_no_tuple_or_case_alias',all(x in disk for x in wanted) and 'RS1' not in disk and (1,43,'A','C') not in disk);disk.close()
after={p:sha(p) for p in BEFORE};check('all_source_runtime_and_prior_review_identity_unchanged',BEFORE==after)
receipt={'status':'PASS_GLOBAL_METHOD_FIXTURES_ONLY' if all(x['pass'] for x in CHECKS) else 'FAILED_GLOBAL_METHOD_FIXTURE_WITNESS','check_count':len(CHECKS),'checks':CHECKS,'case_count':len(CASES),'cases':CASES,'source_bindings_before':BEFORE,'source_bindings_after':after,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'elapsed_seconds':time.monotonic()-START,'full_adapter_admitted':False,'production_body_reads':0,'reference_body_reads':0,'production_mutex_operations':0,'production_workers_or_fits':0,'executable':sys.executable}
with (R/'core_bounded_harmonizer_methods_controls_receipt_v1.json').open('x') as f:json.dump(receipt,f,indent=2);f.write('\n')
print(json.dumps({key:receipt[key] for key in ['status','check_count','case_count','peak_rss_bytes','elapsed_seconds']}))
sys.exit(0 if receipt['status']=='PASS_GLOBAL_METHOD_FIXTURES_ONLY' else 1)
