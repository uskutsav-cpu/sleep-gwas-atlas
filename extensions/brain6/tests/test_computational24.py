"""New v0.3 contracts. All generated inputs below are explicitly artificial."""
import collections
import copy
import itertools
import math
from pathlib import Path
import sqlite3
import subprocess
import sys
import numpy as np
import pytest
from scipy.stats import hypergeom
from brain6.io import ContractError, read_json, write_json, read_tsv, write_tsv, sha256, file_record
from brain6.artifacts import transaction,verify_artifact
from brain6.stats import bh
from brain6.pathway24 import conditional_overlap,run as pathways
from brain6.family24 import joint_placo,local_family
from brain6.compare24 import compare
from brain6.robustness24 import summarize
from brain6.provenance24 import verify_graph
from brain6.scale24 import validate_scale
from brain6.completion24 import draft,audit,seal,validate_plan
from brain6.pairs import freeze_pairs,rank_pairs


def manifest(tmp_path,data,name='spec.json'):
    p=tmp_path/name;write_json(p,data);return p


@pytest.mark.parametrize('N,M,K',[(n,m,k) for n in range(2,9) for m,k in [(1,1),(n//2,n//2),(n-1,2)]])
def test_exact_overlap_bruteforce(N,M,K):
    draws=list(itertools.combinations(range(N),K))
    counts=[sum(x<M for x in d) for d in draws]
    for x in set(counts):
        got=conditional_overlap([(N,M,K,x)])
        assert got['p']==pytest.approx(sum(c>=x for c in counts)/len(counts),abs=1e-12)
        assert got['expected_loci']==pytest.approx(np.mean(counts))
        assert got['variance']==pytest.approx(np.var(counts))


@pytest.mark.parametrize('seed',range(25))
def test_exact_overlap_multistratum(seed):
    rng=np.random.default_rng(seed); strata=[];support=[(0,1.)]
    for _ in range(3):
        n=int(rng.integers(2,9));m=int(rng.integers(0,n+1));k=int(rng.integers(0,n+1))
        lo=max(0,k-n+m);hi=min(m,k);x=int(rng.integers(lo,hi+1));strata.append((n,m,k,x))
        support=[(a+b,p*float(hypergeom.pmf(b,n,m,k))) for a,p in support for b in range(lo,hi+1)]
    got=conditional_overlap(strata);obs=sum(x[3] for x in strata)
    assert got['p']==pytest.approx(sum(p for k,p in support if k>=obs),abs=1e-12)


@pytest.mark.parametrize('strata',[[],[(0,0,0,0)],[(5,6,2,0)],[(5,2,6,0)],[(5,2,2,3)],[(5,5,2,0)],[(5.,2,2,1)]])
def test_exact_overlap_invalid(strata):
    with pytest.raises(ContractError):conditional_overlap(strata)


def units(tmp_path,pairs=('p1','p2'),n=30):
    rs=[]
    for pair in pairs:
        for i in range(n):
            rs.append({'pair_id':pair,'locus_id':f'L{i}','stratum':f's{i%2}','selected':int(i<6),
                'eligible':1,'genes':f'G{i};G{i+100}','genome_build':'GRCh37'})
    p=tmp_path/'units.tsv';write_tsv(p,list(rs[0]),rs);return p,rs


def test_pathway_missing_family_slot(tmp_path):
    p,_=units(tmp_path,pairs=('p1',));gmt=tmp_path/'sets.gmt';gmt.write_text('P1\tsource\tG0\tG1\tG2\nP2\tsource\tG9\tG10\n')
    c={'reviewed':True,'unit':'independent_ld_block','units':str(p),'gmt':str(gmt),'min_genes':1,'max_genes':50,
       'expected_pairs':['p1','p2'],'synthetic':True}
    out=pathways(manifest(tmp_path,c),tmp_path,'a');rs=list(read_tsv(out/'pathways.tsv'))
    assert len(rs)==4
    assert sum(r['status']=='TESTED' for r in rs)==2
    assert read_json(out/'status.json')['family_size']==4
    assert verify_artifact(out)['scientific_status']=='INSUFFICIENT_EVIDENCE'
    got=[float(r['family_fdr']) if r['family_fdr']!='NA' else None for r in rs]
    expected=bh([float(r['p']) if r['p']!='NA' else None for r in rs])
    assert got==pytest.approx(expected,nan_ok=True)


def test_pathway_duplicate_block(tmp_path):
    p,rs=units(tmp_path);write_tsv(p,list(rs[0]),rs+[rs[0]],immutable=False)
    g=tmp_path/'x.gmt';g.write_text('x\tsource\tG1\n')
    c={'reviewed':True,'unit':'independent_ld_block','units':str(p),'gmt':str(g),'expected_pairs':['p1','p2'],'synthetic':True}
    with pytest.raises(ContractError):pathways(manifest(tmp_path,c),tmp_path,'a')


def placo_fixture(tmp_path):
    sources=[];allp=[]
    for pair in ['p1','p2']:
        rs=[]
        for i,p in enumerate([.001,.001,.03,.9,None,None]):
            st='TESTED' if p is not None else 'NUMERICAL_FAILURE' if i==4 else 'EXCLUDED_EXTREME_Z'
            rs.append({'SNP':f'rs{i}','CHR':1,'BP':i+1,'P_PLACO':p if p is not None else 'NA','status':st})
            if st!='EXCLUDED_EXTREME_Z':allp.append(p)
        path=tmp_path/(pair+'.tsv');write_tsv(path,list(rs[0]),rs)
        sources.append({'pair_id':pair,'path':str(path),'sha256':sha256(path),'expected_rows':6,'method':'PLACO_PLUS'})
    c={'synthetic':True,'reviewed':True,'analysis_scope':'GENOMEWIDE_ALL_SELECTED_PAIRS','expected_pairs':['p1','p2'],
       'results':sources,'maximum_failure_rate':.25}
    return c,allp


def test_joint_bh_exact_denominator(tmp_path):
    c,p=placo_fixture(tmp_path);out=joint_placo(manifest(tmp_path,c),tmp_path,'a');rs=list(read_tsv(out/'joint_results.tsv.gz'))
    assert len(rs)==12
    tested=[r for r in rs if r['status']!='EXCLUDED_EXTREME_Z']
    expected=bh(p)
    assert [float(r['JOINT_FAMILY_BH']) if r['JOINT_FAMILY_BH']!='NA' else None for r in tested]==pytest.approx(expected,nan_ok=True)
    status=read_json(out/'status.json');assert status['denominator']==10
    assert status['headline_threshold']==pytest.approx(2.5e-8)


@pytest.mark.parametrize('issue',['missing_pair','short_rows','wrong_hash','zero_p','duplicate_snp','high_failure','mixed_method','candidate_only'])
def test_joint_family_rejects(tmp_path,issue):
    c,_=placo_fixture(tmp_path)
    if issue=='missing_pair':c['results'].pop()
    elif issue=='short_rows':c['results'][0]['expected_rows']=7
    elif issue=='wrong_hash':c['results'][0]['sha256']='f'*64
    elif issue=='candidate_only':c['analysis_scope']='CANDIDATE_ONLY'
    elif issue=='mixed_method':c['results'][0]['method']='PLACO'
    elif issue=='high_failure':c['maximum_failure_rate']=.01
    else:
        src=c['results'][0];p=Path(src['path']);rows=list(read_tsv(p))
        if issue=='zero_p':rows[0]['P_PLACO']=0
        else:rows[1]['SNP']=rows[0]['SNP']
        write_tsv(p,list(rows[0]),rows,immutable=False);src['sha256']=sha256(p)
    f=manifest(tmp_path,c)
    if issue=='high_failure':
        out=joint_placo(f,tmp_path,'a');assert verify_artifact(out)['scientific_status']=='FAILED_QC_NOT_CONSUMED'
    else:
        with pytest.raises((ContractError,sqlite3.IntegrityError)):joint_placo(f,tmp_path,'a')
        assert not (tmp_path/'a').exists()


def test_joint_empirical_needs_native_receipt(tmp_path):
    c,_=placo_fixture(tmp_path);c['synthetic']=False
    with pytest.raises(ContractError):joint_placo(manifest(tmp_path,c),tmp_path,'a')


def test_local_missing_underpowered_distinct(tmp_path):
    p=tmp_path/'plan.tsv';r=tmp_path/'result.tsv'
    write_tsv(p,['pair_id','locus_id'],[dict(pair_id='p',locus_id=str(i)) for i in range(3)])
    write_tsv(r,['pair_id','locus_id','status','p','local_rg','extra'],[
      dict(pair_id='p',locus_id='0',status='TESTED',p=.001,local_rg=.5,extra='keep-safe'),
      dict(pair_id='p',locus_id='1',status='UNIVARIATE_UNDERPOWERED',p='NA',local_rg='NA',extra='ignored')])
    out=local_family(manifest(tmp_path,{'planned_units':str(p),'results':str(r),'reviewed':True,'synthetic':True}),tmp_path,'a')
    s=read_json(out/'status.json');assert s['execution_failures_or_missing']==s['univariate_underpowered']==1
    rs=list(read_tsv(out/'local_family.tsv'));assert float(rs[0]['family_fdr'])==pytest.approx(.003)
    assert rs[-1]['status']=='NOT_RUN'


def comparison(tmp_path,**kw):
    p,_=units(tmp_path);c={'units':str(p),'reviewed':True,'unit':'independent_ld_block',
       'pair_sleep_traits':{'p1':'insomnia','p2':'insomnia'},'synthetic':True,**kw}
    return compare(manifest(tmp_path,c),tmp_path,'a')


def test_shared_sleep_overlap_no_naive_p(tmp_path):
    out=comparison(tmp_path);rs=list(read_tsv(out/'comparisons.tsv'))
    assert rs[0]['shared_blocks']=='6' and rs[0]['p']=='NA'
    assert rs[0]['status']=='DESCRIPTIVE_OVERLAP'


def test_exchangeability_requires_explicit_review(tmp_path):
    with pytest.raises(ContractError):comparison(tmp_path,dependence_strategy='EXCHANGEABILITY_REVIEWED')


def test_reviewed_exchangeability(tmp_path):
    out=comparison(tmp_path,dependence_strategy='EXCHANGEABILITY_REVIEWED',shared_gwas_dependence_justification='ARTIFICIAL independent labels, not a real shared-GWAS claim')
    row=next(read_tsv(out/'comparisons.tsv'));assert float(row['p'])<.01


def test_comparison_missing_pair_preserves_contrasts(tmp_path):
    out=comparison(tmp_path,pair_sleep_traits={'p1':'s','p2':'s','p3':'s'})
    rs=list(read_tsv(out/'comparisons.tsv'));assert len(rs)==3
    assert sum(r['status']=='INSUFFICIENT_COMMON_COVERAGE' for r in rs)==2


@pytest.mark.parametrize('change',['duplicate','build','stratum','unreviewed'])
def test_comparison_invalid(tmp_path,change):
    p,rs=units(tmp_path)
    if change=='duplicate':rs.append(rs[0])
    elif change=='build':rs[-1]['genome_build']='GRCh38'
    elif change=='stratum':rs[-1]['stratum']='wrong'
    write_tsv(p,list(rs[0]),rs,immutable=False)
    c={'units':str(p),'reviewed':change!='unreviewed','unit':'independent_ld_block','pair_sleep_traits':{'p1':'x','p2':'x'},'synthetic':True}
    with pytest.raises(ContractError):compare(manifest(tmp_path,c),tmp_path,'a')


def test_robustness_negative_consistency_not_positive(tmp_path):
    g=tmp_path/'g.tsv';r=tmp_path/'r.tsv'
    write_tsv(g,['unit_id','specification_id','primary','metric','threshold','direction'],[
      dict(unit_id='x',specification_id=str(i),primary=int(i==0),metric='H4',threshold=.8,direction='greater_equal') for i in range(3)])
    write_tsv(r,['unit_id','specification_id','status','value'],[
      dict(unit_id='x',specification_id=str(i),status='PASS',value=.2+i*.01) for i in range(3)])
    c={'reviewed':True,'grid_defined_before_results':True,'grid':str(g),'results':str(r),'synthetic':True}
    out=summarize(manifest(tmp_path,c),tmp_path,'a');rs=list(read_tsv(out/'robustness.tsv'))
    assert rs[0]['status']=='CONSISTENT_WITHIN_TESTED_GRID' and rs[0]['criterion_met_in_all']=='False'


def test_robustness_missing_not_dropped(tmp_path):
    g=tmp_path/'g.tsv';r=tmp_path/'r.tsv'
    write_tsv(g,['unit_id','specification_id','primary','metric','threshold','direction'],[
      dict(unit_id='x',specification_id=str(i),primary=int(i==0),metric='H4',threshold=.8,direction='greater_equal') for i in range(2)])
    write_tsv(r,['unit_id','specification_id','status','value'],[dict(unit_id='x',specification_id='0',status='PASS',value=.9)])
    c={'reviewed':True,'grid_defined_before_results':True,'grid':str(g),'results':str(r),'synthetic':True}
    out=summarize(manifest(tmp_path,c),tmp_path,'a')
    assert next(read_tsv(out/'robustness.tsv'))['missing_or_failed']=='1'


def test_transaction_input_mutation_blocked(tmp_path):
    f=tmp_path/'in';f.write_text('a')
    with pytest.raises(ContractError):
        with transaction(tmp_path,'a',stage='x',inputs=[f],parameters={},synthetic=True) as (w,m):
            (w/'out').write_text('data');f.write_text('b')
    assert not (tmp_path/'a').exists() and any((tmp_path/'_failed').iterdir())


def test_transitive_graph_detects_grandparent_change(tmp_path):
    f=tmp_path/'in';f.write_text('a')
    with transaction(tmp_path,'a',stage='x',inputs=[f],parameters={},synthetic=True) as (w,m):(w/'out').write_text('a')
    with transaction(tmp_path,'b',stage='y',inputs=[tmp_path/'a/receipt.json'],parameters={},synthetic=True) as (w,m):(w/'out').write_text('b')
    roots=[{'path':str(tmp_path/'b'),'fingerprint':verify_artifact(tmp_path/'b')['fingerprint']}]
    assert len(verify_graph(roots,synthetic=True))==2
    f.write_text('b')
    with pytest.raises(ContractError):verify_graph(roots,synthetic=True)


def test_synthetic_graph_cannot_be_empirical(tmp_path):
    with transaction(tmp_path,'a',stage='x',inputs=[],parameters={},synthetic=True) as (w,m):(w/'x').write_text('x')
    with pytest.raises(ContractError):verify_graph([{'path':str(tmp_path/'a'),'fingerprint':None}],synthetic=False)


def locked_plan(tmp_path,atlas_factory):
    m,c,d,*_=atlas_factory();lp=tmp_path/'lock.json';freeze_pairs(m,c,d,lp,'Synthetic reviewer',synthetic=True)
    out=tmp_path/'plan.json';v=draft(lp,out);return out,v


def test_24_missing_scope_blocks_release(tmp_path,atlas_factory):
    f,c=locked_plan(tmp_path,atlas_factory);out=audit(f,tmp_path,'audit',allow_synthetic=True)
    s=read_json(out/'status.json');assert s['unresolved']>100 and not s['publication_ready']
    rows=list(read_tsv(out/'steps_1_24.tsv'));assert len(rows)==24 and rows[-1]['status']=='BLOCKED'
    with pytest.raises(ContractError):seal(f,out,tmp_path,'release',allow_synthetic=True)


@pytest.mark.parametrize('issue',['missing_step','disable_locus_expansion','wrong_stage','wrong_method','unreviewed_omission','changed_matrix'])
def test_24_plan_integrity(tmp_path,atlas_factory,issue):
    f,c=locked_plan(tmp_path,atlas_factory)
    if issue=='missing_step':c['requirements'].pop()
    elif issue=='disable_locus_expansion':next(r for r in c['requirements'] if r['step']==10)['expanded_scope_required']=False
    elif issue=='wrong_stage':next(r for r in c['requirements'] if r['step']==17)['accepted_artifact_stages']=['normalize']
    elif issue=='wrong_method':next(r for r in c['requirements'] if r['step']==17)['accepted_methods']=['h2']
    elif issue=='unreviewed_omission':c['requirements'][0]['required']=False
    elif issue=='changed_matrix':(tmp_path/'matrix.tsv').write_text('changed')
    with pytest.raises(ContractError):validate_plan(c)


def test_24_one_artifact_not_all_loci(tmp_path,atlas_factory):
    f,c=locked_plan(tmp_path,atlas_factory);scope=tmp_path/'scope.tsv';write_tsv(scope,['unit_id'],[{'unit_id':'l1'},{'unit_id':'l2'}])
    with transaction(tmp_path,'ld',stage='native_job',inputs=[],parameters={'method':'ld'},synthetic=True) as (w,m):
        (w/'x').write_text('x');m['scientific_status']='PASS'
    req=next(r for r in c['requirements'] if r['step']==10);req['scope_manifest']=file_record(scope)
    rr=verify_artifact(tmp_path/'ld');req['bindings']=[{'unit_id':u,'path':str(tmp_path/'ld'),'fingerprint':rr['fingerprint']} for u in ['l1','l2']]
    write_json(f,c,immutable=False)
    with pytest.raises(ContractError):audit(f,tmp_path,'audit',allow_synthetic=True)


def test_scale_binary_lmm_is_not_logodds():
    c=dict(scale_reviewed=True,n_semantics='effective',sample_size_justification='source',trait_type='cc',effect_scale='binary_lmm_beta')
    with pytest.raises(ContractError):validate_scale('susie',c)
    c.update(trait_type='quant',sdY_source='Source sample variance')
    with pytest.raises(ContractError):validate_scale('susie',c)
    c.update(binary_lmm_approximation_reviewed=True,binary_lmm_approximation_justification='Explicit model approximation')
    validate_scale('susie',c)


def test_scale_mr_prevents_fake_oddsratio():
    c=dict(scale_reviewed=True,exposure_effect_scale='log_odds',exposure_effect_unit='logodds',outcome_effect_scale='quantitative_beta',outcome_effect_unit='hours',report_odds_ratios=True)
    with pytest.raises(ContractError):validate_scale('mr',c)
    c['report_odds_ratios']=False;validate_scale('mr',c)


@pytest.mark.parametrize('method,cfg',[('mr',{}),('susie',{}),('susie',dict(scale_reviewed=True,trait_type='quant',effect_scale='quantitative_beta'))])
def test_scales_no_unsourced_defaults(method,cfg):
    with pytest.raises(ContractError):validate_scale(method,cfg)


def test_rank_nan_cannot_be_selected(atlas_factory):
    m,c,d,rows,cfg,dec=atlas_factory();rows[0]['fdr']='nan';write_tsv(m,list(rows[0]),rows,immutable=False)
    target=next(x for x in rank_pairs(m,cfg) if x['pair_id']=='insomnia__adhd')
    assert target['fdr'] is None and not target['eligible']


def test_cell_family_keeps_missing_test_slots(tmp_path):
    from brain6.cell24 import run
    score=tmp_path/'scores.tsv'
    rows=[dict(pair_id='p1',cell_type='cell1',locus_id=f'L{i}',score=1 if i<5 else 0,selected=int(i<5),stratum='s') for i in range(20)]
    write_tsv(score,list(rows[0]),rows)
    c=dict(reviewed=True,unit='independent_ld_block',expected_pairs=['p1','p2'],expected_cells=['cell1','cell2'],
        scores=str(score),synthetic=True,permutations=99,seed=92,exchangeability_justification='Synthetic matched blocks')
    out=run(manifest(tmp_path,c),tmp_path,'a');r=list(read_tsv(out/'cell_enrichment.tsv'))
    assert len(r)==4 and sum(x['status']=='MISSING_SCORES' for x in r)==3
    assert float(r[0]['family_fdr'])==pytest.approx(min(1,float(r[0]['p'])*4))
    assert verify_artifact(out)['scientific_status']=='INSUFFICIENT_EVIDENCE'


def test_validation_executes_real_test_and_preserves_skip(tmp_path):
    from brain6.validation24 import run
    test=tmp_path/'test_artificial_example.py'
    test.write_text('import pytest\ndef test_true(): assert 1+1==2\ndef test_skipped(): pytest.skip("explicit unavailable test")\n')
    c=dict(reviewed=True,synthetic=True,mode='regression',test_suites=[dict(id='demo',path=str(test))],timeout_seconds=30)
    out=run(manifest(tmp_path,c),tmp_path,'a');status=read_json(out/'status.json')
    assert status['passed']==status['skipped']==1
    assert status['status']=='INSUFFICIENT_EVIDENCE'


def test_native_validation_no_arbitrary_fake_test(tmp_path):
    from brain6.validation24 import run
    f=tmp_path/'test_true.py';f.write_text('def test_x(): assert True\n')
    c=dict(reviewed=True,synthetic=True,mode='native_environment',test_suites=[dict(id='fake',path=str(f))])
    with pytest.raises(ContractError):run(manifest(tmp_path,c),tmp_path,'a')


def test_junit_cannot_be_empty(tmp_path):
    from brain6.validation24 import junit_counts
    p=tmp_path/'empty.xml';p.write_text('<testsuites/>')
    with pytest.raises(ContractError):junit_counts(p)



def test_24_covariance_is_not_a_validated_model(tmp_path,atlas_factory):
    f,c=locked_plan(tmp_path,atlas_factory)
    settings=tmp_path/'cov.settings.json';write_json(settings,{'mode':'covariance'})
    params={'method':'genomicsem','inputs':{'settings':file_record(settings)}}
    with transaction(tmp_path,'cov',stage='native_job',inputs=[settings],parameters=params,synthetic=True) as (w,m):
        (w/'S.tsv').write_text('ARTIFICIAL covariance fixture');m['scientific_status']='PASS'
    req=next(r for r in c['requirements'] if r['step']==18)
    req['bindings']=[{'unit_id':'WHOLE_SCOPE','path':str(tmp_path/'cov'),'fingerprint':verify_artifact(tmp_path/'cov')['fingerprint']}]
    write_json(f,c,immutable=False)
    with pytest.raises(ContractError,match='Covariance preparation'):
        audit(f,tmp_path,'audit',allow_synthetic=True)
