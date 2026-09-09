from pathlib import Path
import pytest
from brain6.io import write_json,write_tsv,read_json,read_tsv,sha256,ContractError
from brain6.artifacts import transaction,verify_artifact,verify_current_artifact
from brain6.pairs import freeze_pairs
from brain6.pipeline import run_pipeline
from brain6.gwas import normalize,join_pair
from brain6.ld import index_pair
from brain6.pleiotropy import collate_placo
from brain6.followup import run_followup


def test_preparation_dag_executes_and_resumes(atlas_factory,source_factory,tmp_path):
    m,c,d,*_=atlas_factory();lock=tmp_path/'lock.json';freeze_pairs(m,c,d,lock,'tester',synthetic=True)
    a,_,_=source_factory();b,_,_=source_factory()
    runtime={'schema_version':1,'reviewed':True,'synthetic':True,'pair_lock':str(lock),
       'output_root':str(tmp_path/'run'),'tasks':[
        {'id':'a','kind':'normalize','parameters':{'source':str(a)}},
        {'id':'b','kind':'normalize','parameters':{'source':str(b)}},
        {'id':'pair','kind':'join','depends_on':['a','b'],'parameters':{'left':'@a','right':'@b','min_variants':1}},
        {'id':'index','kind':'index','depends_on':['pair'],'parameters':{'pair_dir':'@pair'}}]}
    path=tmp_path/'runtime.json';write_json(path,runtime)
    run_pipeline(path,allow_synthetic=True);r=verify_artifact(tmp_path/'run/index')
    run_pipeline(path,allow_synthetic=True);assert r==verify_artifact(tmp_path/'run/index')
    runtime['tasks'][2]['parameters']['min_variants']=2;write_json(path,runtime,immutable=False)
    with pytest.raises(ContractError,match='changed'):run_pipeline(path,allow_synthetic=True)


def test_current_artifact_detects_input_change(tmp_path):
    p=tmp_path/'input';p.write_text('original')
    with transaction(tmp_path,'artifact',stage='test',inputs=[p],parameters={},synthetic=True) as (w,m):
        (w/'output').write_text('test')
    verify_current_artifact(tmp_path/'artifact');p.write_text('changed')
    with pytest.raises(ContractError,match='Stale'):verify_current_artifact(tmp_path/'artifact')


def test_failed_transaction_preserves_diagnostics(tmp_path):
    with pytest.raises(ValueError):
        with transaction(tmp_path,'bad',stage='test',inputs=[],parameters={},synthetic=True) as (w,m):
            (w/'native.log').write_text('a useful failure diagnostic');raise ValueError('test failure')
    assert not (tmp_path/'bad').exists()
    failed=next((tmp_path/'_failed').iterdir())
    assert (failed/'native.log').read_text()=='a useful failure diagnostic'
    assert read_json(failed/'receipt.json')['status']=='FAILED'

@pytest.mark.parametrize('has_locus',[False,True])
def test_followup_accounts_for_zero_hits_and_blocked_ld(source_factory,atlas_factory,tmp_path,monkeypatch,has_locus):
    # Every native response below is an explicitly synthetic test double, NOT genetics execution.
    a,_,_=source_factory();b,_,_=source_factory()
    l=normalize(a,tmp_path,'l',synthetic=True);r=normalize(b,tmp_path,'r',synthetic=True)
    pair=join_pair(l,r,tmp_path,'pair',min_variants=1,synthetic=True)
    index=index_pair(pair,tmp_path,'index')
    rows=[dict(SNP='rs1',CHR=1,BP=100,Z1=2,Z2=3,P_PLACO=1e-10,status='TESTED')]
    f=tmp_path/'fixture_chunk.tsv';write_tsv(f,list(rows[0]),rows)
    manifest=tmp_path/'chunks.tsv';write_tsv(manifest,['path','sha256'],[{'path':str(f),'sha256':sha256(f)}])
    col=collate_placo(manifest,tmp_path,'collated',expected_rows=1,n_pairs=6,synthetic=True)
    m,c,d,*_=atlas_factory();lock=tmp_path/'lock.json';freeze_pairs(m,c,d,lock,'tester',synthetic=True)
    prefix=str(tmp_path/'reference')
    for s in ['.bed','.bim','.fam']:Path(prefix+s).write_text('synthetic fake reference, not genotype data')
    blocks=tmp_path/'blocks.tsv';write_tsv(blocks,['LOC','CHR','START','STOP'],[dict(LOC='b1',CHR=1,START=1,STOP=1000)])
    settings={'reviewed':True,'method_compatibility_reviewed':True,'pair_id':'insomnia__mdd',
       'pair_lock':str(lock),'pair_index':str(index),'collated':str(col),
       'discovery_method':'PLACO_PLUS','ancestry':'EUR','genome_build':'GRCh37',
       'maximum_locus_failure_rate':.01,'reference_prefix':prefix,'ld_blocks':str(blocks),
       'output_root':str(tmp_path/'deep'),'clumping':{'plink':'plink','r2':.1,'kb':500},
       'max_working_bytes':10**8,'minimum_ld_coverage':.8}
    sp=tmp_path/'deep.json';write_json(sp,settings)
    def fake_native(job,root,**kwargs):
        j=read_json(job);name=j['job_id']
        if name!='clump':raise ContractError('SYNTHETIC BLOCKED_BY_SOFTWARE: no native LD engine')
        target=Path(root)/name
        if target.exists():return verify_artifact(target)
        with transaction(root,name,stage='SYNTHETIC_NATIVE_STUB',inputs=[],parameters={},synthetic=True) as (w,meta):
            leads=[dict(SNP='rs1',CHR=1,BP=100,P=1e-10)] if has_locus else []
            write_tsv(w/'leads.tsv',['SNP','CHR','BP','P'],leads)
            meta['scientific_status']='PASS' if has_locus else 'NO_SIGNAL'
        return verify_artifact(target)
    monkeypatch.setattr('brain6.followup.run_job',fake_native)
    result=run_followup(sp,allow_synthetic=True)
    assert result['scientific_status']==('FAILED_QC_NOT_CONSUMED' if has_locus else 'NO_SIGNAL')
    summary=read_json(tmp_path/'deep/summary/summary.json')
    assert summary['failures']==int(has_locus)
    assert run_followup(sp,allow_synthetic=True)['fingerprint']==result['fingerprint']
