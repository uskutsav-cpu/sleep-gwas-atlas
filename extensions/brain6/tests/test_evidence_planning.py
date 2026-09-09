import json
import sqlite3
from pathlib import Path
import pytest
from brain6.io import ContractError,write_tsv,write_json,read_json,read_tsv,sha256
from brain6.evidence import integrate,cell_enrichment,cross_disorder
from brain6.pairs import freeze_pairs,write_ranking,BRAIN
from brain6.assembly import assemble_preparation
from brain6.bindings import discover_inputs
from brain6.matlab import matlab_string,config_entries
from brain6.planning import placo_jobs,execute_placo_plan
from brain6.gwas import normalize,join_pair
from brain6.pleiotropy import split_pair
from brain6.pipeline import create_job


def make_evidence(tmp_path,etype='nearest_gene',h4='.9',molecular_h4='NA',status='PASS'):
    a=tmp_path/'loci.tsv';b=tmp_path/'mol.tsv'
    write_tsv(a,['pair_id','locus_id','status','trait_h4','pip_max'],[dict(pair_id='insomnia__mdd',locus_id='b1',status=status,trait_h4=h4,pip_max=.9)])
    row=dict(pair_id='insomnia__mdd',locus_id='b1',gene_id='GENE',evidence_type=etype,molecular_h4=molecular_h4,tissue='brain',cell_type='cell',source='fixture_source')
    write_tsv(b,list(row),[row]);return a,b

@pytest.mark.parametrize('etype,h4,mh,status,expected',[
 ('nearest_gene','.9','NA','PASS','SHARED_SIGNAL_GENE_UNRESOLVED'),
 ('expression','.9','NA','PASS','SHARED_SIGNAL_GENE_UNRESOLVED'),
 ('eqtl_coloc','.9','.95','PASS','SHARED_SIGNAL_WITH_MOLECULAR_SUPPORT'),
 ('coding','.9','NA','PASS','SHARED_SIGNAL_WITH_CODING_ANNOTATION'),
 ('eqtl_coloc','.2','.95','PASS','REGIONAL_OR_DISTINCT_SIGNAL'),
 ('eqtl_coloc','.9','.95','FAILED','UNRESOLVED_OR_BLOCKED')])
def test_evidence_tiers(tmp_path,etype,h4,mh,status,expected):
    a,b=make_evidence(tmp_path,etype,h4,mh,status)
    rows=integrate(a,b,tmp_path/'out.tsv');assert rows[0]['tier']==expected

@pytest.mark.parametrize('h4',['-1','nan','inf','1.2'])
def test_bad_evidence_probabilities(tmp_path,h4):
    a,b=make_evidence(tmp_path,h4=h4)
    with pytest.raises(ContractError):integrate(a,b,tmp_path/'o.tsv')

def test_empty_molecular_preserves_locus(tmp_path):
    a,b=make_evidence(tmp_path)
    write_tsv(b,['pair_id','locus_id','gene_id','evidence_type','molecular_h4','tissue','cell_type','source'],[],immutable=False)
    assert len(integrate(a,b,tmp_path/'o.tsv'))==1

def test_locus_enrichment_reproducible(tmp_path):
    rows=[dict(pair_id='p',cell_type=c,locus_id=f'l{i}',score=i,selected=int(i>=3),stratum='one') for c in ['a','b'] for i in range(6)]
    p=tmp_path/'scores.tsv';write_tsv(p,list(rows[0]),rows)
    a=cell_enrichment(p,tmp_path/'a.tsv',permutations=99);b=cell_enrichment(p,tmp_path/'b.tsv',permutations=99)
    assert a==b and len(a)==2 and all(0<=r['family_fdr']<=1 for r in a)

def test_reject_double_counted_locus(tmp_path):
    row=dict(pair_id='p',cell_type='c',locus_id='l',score=2,selected=1,stratum='s')
    p=tmp_path/'s.tsv';write_tsv(p,list(row),[row,row])
    with pytest.raises(ContractError):cell_enrichment(p,tmp_path/'o.tsv',permutations=99)

def test_reference_prefix_checks_all_members(tmp_path):
    prefix=str(tmp_path/'ref')
    for suffix in ['.bed','.bim','.fam']:Path(prefix+suffix).write_text('fixture')
    d=discover_inputs('ld',{'reference_prefix':prefix});assert len(d)==3
    Path(prefix+'.fam').unlink()
    with pytest.raises(ContractError):discover_inputs('ld',{'reference_prefix':prefix})

def test_reference_prefix_no_guess(tmp_path):
    with pytest.raises(ContractError):discover_inputs('ld',{'reference_prefix':'UNRESOLVED'})

def test_r_job_binds_inputs_and_source(tmp_path):
    p=tmp_path/'fit.rds';p.write_bytes(b'not real R data')
    s=tmp_path/'s.json';write_json(s,{'fit1':str(p),'fit2':str(p)})
    job=create_job(s,method='coloc',job_id='one',inputs={},outputs={'semantic_status':'status.json'},root=tmp_path/'jobs',reviewed=True,synthetic=True)
    j=read_json(job);assert {'adapter','settings','common'}<=set(j['inputs'])
    assert any(r['sha256']==sha256(p) for r in j['inputs'].values())
    assert j['argv'][0]=='Rscript'

def test_matlab_quote_and_config(tmp_path):
    assert matlab_string("a'b")=="'a''b'"
    with pytest.raises(ContractError):matlab_string('a\nb')
    p=tmp_path/'c.txt';p.write_text('# comment\nstat=one\nx=2\n')
    assert config_entries(p)=={'stat':'one','x':'2'}
    p.write_text('x=1\nx=2\n')
    with pytest.raises(ContractError):config_entries(p)

def test_assembly_deduplicates_source_tasks(atlas_factory,tmp_path):
    m,c,d,rows,cfg,decs=atlas_factory();lp=tmp_path/'lock.json';freeze_pairs(m,c,d,lp,'tester',synthetic=True)
    sources={}
    for t in ['insomnia']+BRAIN:
        s=tmp_path/(t+'.json');write_json(s,{'trait_id':t,'synthetic':True});sources[t]=str(s)
    sf=tmp_path/'sources.json';write_json(sf,sources)
    runtime=assemble_preparation(lp,sf,tmp_path/'work',tmp_path/'runtime.json',reviewed=True)
    assert sum(t['kind']=='normalize' for t in runtime['tasks'])==7
    assert sum(t['kind']=='join' for t in runtime['tasks'])==6
    assert not any(t['kind']=='split' and t['pair_id']=='insomnia__adhd' for t in runtime['tasks'])

def test_assembly_missing_sources(atlas_factory,tmp_path):
    m,c,d,*_=atlas_factory();lp=tmp_path/'lock.json';freeze_pairs(m,c,d,lp,'tester',synthetic=True)
    s=tmp_path/'s.json';write_json(s,{})
    with pytest.raises(ContractError):assemble_preparation(lp,s,tmp_path/'w',tmp_path/'r',reviewed=True)

def test_placo_planner_protects_adhd(tmp_path):
    with pytest.raises(ContractError,match='read-only'):
        placo_jobs('missing','missing','missing','insomnia__adhd','missing','missing',tmp_path/'plan',reviewed=True)

def test_placo_planner_real_pair_chunks(source_factory,atlas_factory,tmp_path):
    a,_,_=source_factory();b,_,_=source_factory()
    left=normalize(a,tmp_path,'left',synthetic=True);right=normalize(b,tmp_path,'right',synthetic=True)
    pair=join_pair(left,right,tmp_path,'pair',min_variants=1,synthetic=True)
    chunks=split_pair(pair,tmp_path,'chunks',chunk_rows=3)
    m,c,d,*_=atlas_factory();lp=tmp_path/'pairlock.json';freeze_pairs(m,c,d,lp,'tester',synthetic=True)
    src=tmp_path/'official.R';src.write_text('# SYNTHETIC source, not PLACO')
    plan=placo_jobs(pair,chunks,lp,'insomnia__mdd',src,sha256(src),tmp_path/'plan',reviewed=True)
    assert read_json(plan)['n_pairs']==6
    data=read_json(plan);data['n_pairs']=1;write_json(plan,data,immutable=False)
    with pytest.raises(ContractError,match='changed'):execute_placo_plan(plan,allow_synthetic=True)

def test_figures_use_full_matrix_and_actual_ps(atlas_factory,tmp_path):
    from brain6.plots import heatmap,manhattan,qq
    m,c,*_=atlas_factory();r=tmp_path/'SYNTHETIC_ranking.tsv';write_ranking(m,c,r)
    heatmap(r,tmp_path/'heat.png');assert (tmp_path/'heat.png').stat().st_size>1000
    values=[dict(CHR=1,BP=i*10+1,P_PLACO=.001*(i+1),status='TESTED') for i in range(10)]
    p=tmp_path/'p.tsv';write_tsv(p,list(values[0]),values)
    result=manhattan(p,tmp_path/'man.png');assert result['input']==10
    db=sqlite3.connect(tmp_path/'r.sqlite');db.execute('CREATE TABLE results(p REAL,status TEXT)');db.executemany('INSERT INTO results VALUES(?,?)',[(r['P_PLACO'],r['status']) for r in values]);db.commit();db.close()
    qq(tmp_path/'r.sqlite',tmp_path/'qq.png');assert (tmp_path/'qq.png').stat().st_size>1000
