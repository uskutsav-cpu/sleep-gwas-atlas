#!/usr/bin/env python3
"""Independent synthetic-only identity/math oracle; no body or asset reads."""
from decimal import Decimal, getcontext
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys

PACKAGE=Path(__file__).resolve().parents[1]
ROOT=PACKAGE.parent
SCRIPTS=PACKAGE/'scripts'
sys.path.insert(0,str(SCRIPTS))
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for data in iter(lambda:f.read(65536),b''):h.update(data)
    return h.hexdigest()
def module(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    value=importlib.util.module_from_spec(spec);spec.loader.exec_module(value);return value
paths=[SCRIPTS/'56_preprocess_finngen_insomnia_feasibility.py',SCRIPTS/'56_preprocess_finngen_insomnia_feasibility_v2.py',SCRIPTS/'57_verify_finngen_row_controls.py',SCRIPTS/'57_verify_finngen_row_controls_v2.py',ROOT/'scripts/liftover_chain.py',ROOT/'discovery_extension/scripts/47_stream_replication_sources.py',PACKAGE/'FROZEN_FINNGEN_INSOMNIA_FEASIBILITY_PROTOCOL_v1.md',PACKAGE/'FROZEN_FINNGEN_PREPROCESSING_DETAIL_v1.md',PACKAGE/'logs/finngen_row_controls_v4.json',PACKAGE/'logs/finngen_row_controls_v4_2.json']
before={str(p):sha(p) for p in paths}
old=module(paths[0],'original_finngen_row_only')
new=module(paths[1],'corrected_finngen_row_only')
lift=module(ROOT/'scripts/liftover_chain.py','independent_finngen_synthetic_chain')
chain=lift.ChainIndex({1:[(100,110,1,200,1000,'+',1),(120,130,1,300,1000,'-',2)]})
lookup={'rs1':('A','C',1,201),'rs2':('T','G',1,700)}
base={'#chrom':'1','pos':'101','ref':'C','alt':'A','rsids':'rs1','nearest_genes':'IGNORED','pval':'0.9','mlogp':'0','beta':'0.2','sebeta':'0.1','af_alt':'0.3','af_alt_cases':'NA','af_alt_controls':'NA'}
checks=[]
def compare(label,actual,expected):
    ok=actual==expected
    checks.append({'label':label,'actual':actual,'expected':expected,'pass':ok})
    if not ok:raise AssertionError(label)
def row(changes=None,mapper=chain,seen=None):return new.qualify_row({**base,**(changes or {})},lookup,mapper,seen or set())
compare('original_upper_boundary_failure_preserved',old.qualify_row({**base,'af_alt':'0.99'},lookup,chain,set())[0],'retained')
for frequency,expected in [('0.01','maf_at_or_below_0_01_or_invalid'),('0.99','maf_at_or_below_0_01_or_invalid'),('0.989','retained'),('0.0101','retained'),('0','maf_at_or_below_0_01_or_invalid'),('1','maf_at_or_below_0_01_or_invalid'),('inf','maf_at_or_below_0_01_or_invalid')]:compare('strict_frequency_'+frequency,row({'af_alt':frequency})[0],expected)
for effect,other,sign in [('A','C',1),('C','A',-1),('T','G',1),('G','T',-1)]:compare('allele_sign_'+effect+other,row({'alt':effect,'ref':other})[1]['Z'],sign*2.0)
compare('reverse_coordinate_1_based',chain.map_point('chr1','121'),('mapped',(1,700,'-')))
compare('reverse_native_beta_not_negated_by_strand',row({'rsids':'rs2','pos':'121'})[1]['Z'],2.0)
compare('reverse_effect_swap_negates_beta',row({'rsids':'rs2','pos':'121','alt':'C','ref':'A'})[1]['Z'],-2.0)
compare('forward_exclusive_end',chain.map_point('1','111'),('unmapped',None))
compare('wrong_coordinate',row({'pos':'102'})[0],'mapped_coordinate_rsid_reference_mismatch')
compare('multiple_distinct_HM3_ids',row({'rsids':'rs1;rs2'})[0],'multiple_hm3_rsids')
compare('same_alias_identity',row({'rsids':'rs1;rs1'})[0],'retained')
compare('retained_only_duplicate',row(seen={'rs1'})[0],'duplicate_retained_hm3_rsid')
for p in ('NA','0','1e-300','2','-1'):compare('P_never_filters_'+p,row({'pval':p})[0],'retained')
for beta,se in [('nan','0.1'),('0.2','0'),('0.2','-1'),('0.2','inf')]:compare('supplied_stat_failure_'+beta+'_'+se,row({'beta':beta,'sebeta':se})[0],'invalid_beta_or_supplied_se')
getcontext().prec=60
exact=Decimal(4)*Decimal(51643)*Decimal(446273)/Decimal(497916)
native_n=4/(1/51643+1/446273)
relative=abs(Decimal.from_float(native_n)-exact)/exact
assert relative<Decimal('1e-15')
controls=[]
for path,count in [(PACKAGE/'logs/finngen_row_controls_v4.json',32),(PACKAGE/'logs/finngen_row_controls_v4_2.json',34)]:
    r=json.loads(path.read_text());compare('bound_control_count_'+path.name,len(r['checks']),count)
    assert r['all_pass'] and all(c['passed'] for c in r['checks'])
    for name,digest in r['dependencies_sha256'].items():assert sha(name)==digest
    controls.append({'path':str(path),'sha256':sha(path),'synthetic_controls':count})
assert {str(p):sha(p) for p in paths}==before
record={'status':'PASS_CORRECTED_SCIENTIFIC_ROW_RULES_NO_SOURCE_ADMISSION',
        'checker_sha256':sha(__file__),'consumed_small_file_sha256':before,'all_checks_pass':True,'checks':checks,
        'initial_strict_MAF_upper_boundary_defect_preserved':True,'corrected_direct_boundaries_pass':True,
        'assumed_constant_effective_N_float':native_n,'exact_decimal_N':str(exact),'float_relative_difference':str(relative),'historical_12g_N':format(native_n,'.12g'),
        'supplied_normal_P_for_Z2':math.erfc(2/math.sqrt(2)),'bound_prior_controls':controls,
        'real_GWAS_or_reference_or_chain_asset_body_reads':0,'estimator_or_real_preprocessing_launches':0,
        'full_stream_source_gzip_and_row_claims_independently_verified':False,'INFO_or_per_variant_N_or_replication_admitted':False}
out=PACKAGE/'reviews/independent_finngen_row_science_receipt_v4.json'
with out.open('x') as f:json.dump(record,f,indent=2);f.write('\n')
print(json.dumps({'status':record['status'],'checks':len(checks),'assumed_N':native_n,'original_bug_preserved':True}))
