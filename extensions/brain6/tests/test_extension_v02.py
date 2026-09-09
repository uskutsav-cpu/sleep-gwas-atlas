"""Real Python tests on explicitly synthetic fixtures; native receipts here are stubs."""
import json
from pathlib import Path
import pytest
from brain6.io import ContractError,read_json,write_json,write_tsv,read_tsv,sha256,file_record
from brain6.artifacts import transaction,verify_artifact,fingerprint
from brain6.pairs import BRAIN,freeze_pairs,verify_lock
from brain6.gwas import normalize,FIELDS
from brain6.mr_workflow import candidates,harmonize,collate
from brain6.replication import evaluate
from brain6.campaign import template,audit,STAGES
from brain6.qtl import index,locus,REQUIRED
from brain6.annotations import IntervalIndex,score
from brain6.pleiotropy import split_pair

@pytest.fixture
def pair_lock(atlas_factory,tmp_path):
    m,c,d,*_=atlas_factory();p=tmp_path/'pairs.json'
    freeze_pairs(m,c,d,p,'Synthetic test fixture',synthetic=True)
    return p

@pytest.mark.parametrize('text',['{"a":1,"a":2}','{"x":NaN}','{"x":Infinity}','{"x":-Infinity}'])
def test_strict_json(tmp_path,text):
    p=tmp_path/'bad.json';p.write_text(text)
    with pytest.raises(ContractError):read_json(p)


def test_all_null_panel_records_all_six(atlas_factory,tmp_path):
    m,c,d,rows,cfg,decisions=atlas_factory()
    for r in rows:r['p']=r['fdr']=.8
    m.unlink();write_tsv(m,list(rows[0]),rows)
    decisions=[dict(disease_trait=b,sleep_trait='',role='NO_ELIGIBLE_PAIR',reason='No eligible synthetic correlation') for b in BRAIN]
    d.unlink();write_tsv(d,list(decisions[0]),decisions)
    p=tmp_path/'null.json';r=freeze_pairs(m,c,d,p,'Synthetic reviewer',synthetic=True)
    assert r['n_pairs']==0 and r['genomewide_family_threshold'] is None
    assert len(verify_lock(p)['disorder_decisions'])==6
    t=tmp_path/'campaign.json';template(p,t)
    out=audit(t,tmp_path/'out','coverage')
    rs=list(read_tsv(out/'coverage.tsv'))
    assert sum(r['status']=='NO_ELIGIBLE_PAIR' for r in rs)==6


def test_native_code_changes_fingerprint(tmp_path,monkeypatch):
    import brain6.artifacts as a
    code=tmp_path/'brain6';(code/'resources').mkdir(parents=True)
    (code/'artifacts.py').write_text('# fixture')
    r=code/'resources/test.R';r.write_text('x <- 1')
    monkeypatch.setattr(a,'__file__',str(code/'artifacts.py'))
    key1,_=fingerprint([],{},'x');r.write_text('x <- 2');key2,_=fingerprint([],{},'x')
    assert key1!=key2


def test_receipt_requires_nonempty_inventory(tmp_path):
    p=tmp_path/'fake';p.mkdir();write_json(p/'receipt.json',{'status':'COMPLETE','outputs':[]})
    with pytest.raises(ContractError,match='inventory'):verify_artifact(p)


def test_mr_exposure_only_thresholds(source_factory,tmp_path):
    spec,_,rows=source_factory()
    rows[0].update(P=1e-12,BETA=.2,SE=.02)
    rows[1].update(P=1e-12,BETA=.01,SE=.02) # large P evidence but weak F: excluded
    s,_,_=source_factory(rows)
    x=normalize(s,tmp_path/'work','x',synthetic=True)
    c=candidates(x,tmp_path/'work','c')
    kept=list(read_tsv(c/'candidates.tsv'))
    assert [r['SNP'] for r in kept]==['rs1']
    assert read_json(c/'status.json')['outcome_used_for_selection'] is False


def test_mr_no_candidates_is_not_failure(source_factory,tmp_path):
    s,_,_=source_factory();x=normalize(s,tmp_path/'w','x',synthetic=True)
    c=candidates(x,tmp_path/'w','c')
    assert verify_artifact(c)['scientific_status']=='NO_SIGNAL'
    assert list(read_tsv(c/'candidates.tsv'))==[]


def synthetic_clump(tmp_path,cand,snps,*,wrong_source=False):
    # This is deliberately NOT a PLINK execution. It tests the downstream contract.
    settings=tmp_path/'clump-settings.json'
    write_json(settings,{'associations':str(cand/'candidates.tsv') if not wrong_source else '/wrong/placo.tsv',
                        'p_column':'P','p1':5e-8,'r2':.001,'kb':10000})
    with transaction(tmp_path/'w','clump',stage='native_job',inputs=[settings],
        parameters={'method':'clump','inputs':{'settings':file_record(settings)}},synthetic=True) as (work,meta):
        write_tsv(work/'leads.tsv',['SNP'],[{'SNP':s} for s in snps]);meta['scientific_status']='PASS'
    return tmp_path/'w/clump'


def mr_fixture(source_factory,tmp_path,wrong_source=False):
    _,_,rows=source_factory()
    for r in rows:r.update(P=1e-10,BETA=.2,SE=.02)
    sx,_,_=source_factory(rows);x=normalize(sx,tmp_path/'w','x',synthetic=True)
    ys=[{**r,'P':.9,'BETA':-.03,'A1':'C','A2':'A','EAF':.7} for r in rows[:8]]
    sy,_,_=source_factory(ys);y=normalize(sy,tmp_path/'w','y',synthetic=True)
    cand=candidates(x,tmp_path/'w','cand')
    clump=synthetic_clump(tmp_path,cand,['rs1','rs2','rs3','rs10'],wrong_source=wrong_source)
    return cand,clump,y


def test_mr_keeps_nonsignificant_outcomes_and_aligns(source_factory,tmp_path):
    c,l,y=mr_fixture(source_factory,tmp_path)
    result=harmonize(c,l,y,tmp_path/'w','mr')
    ys=list(read_tsv(result/'outcome.tsv'));a=list(read_tsv(result/'attrition.tsv'))
    assert len(ys)==3 and all(float(r['P'])==.9 for r in ys)
    assert all(float(r['BETA'])==.03 and r['A1']=='A' for r in ys)
    assert a[-1]['status']=='MISSING_OUTCOME'


def test_mr_rejects_outcome_or_placo_instrument_selection(source_factory,tmp_path):
    c,l,y=mr_fixture(source_factory,tmp_path,True)
    with pytest.raises(ContractError,match='exposure-only'):harmonize(c,l,y,tmp_path/'w','mr')


def test_mr_family_keeps_twelve_slots(pair_lock,tmp_path):
    p=tmp_path/'mr.json';write_json(p,{'pair_lock':str(pair_lock),'reviewed':True,'results':[]})
    result=collate(p,tmp_path/'w','mr')
    rs=list(read_tsv(result/'mr_family.tsv'))
    assert len(rs)==12 and all(r['family_fdr']=='NA' for r in rs)
    assert read_json(result/'status.json')['tested']==0


def test_mr_family_rejects_duplicate_units(pair_lock,tmp_path):
    b={'pair_id':'insomnia__adhd','direction':'FORWARD','status':'UNDERPOWERED','reason':'Synthetic'}
    p=tmp_path/'mr.json';write_json(p,{'pair_lock':str(pair_lock),'reviewed':True,'results':[b,b]})
    with pytest.raises(ContractError,match='duplicate'):collate(p,tmp_path/'w','mr')


def replication_manifest(pair_lock,tmp_path,**changes):
    cohort=tmp_path/'cohorts.md';cohort.write_text('SYNTHETIC reviewed cohort comparison')
    r={'pair_id':'insomnia__adhd','source':'synthetic://replication','cohort_audit':str(cohort),
       'cohort_audit_sha256':sha256(cohort),'independence_scope':'DISEASE_ONLY','phenotype_match':'exact',
       'independence_reviewed':True,'qc_pass':True,'rg':.4,'se':.05,'p':1e-10,
       'sleep_h2_z':8,'disease_h2_z':8,'discovery_cohorts':['D'],'replication_cohorts':['R']}
    r.update(changes)
    p=tmp_path/'rep.json';write_json(p,{'pair_lock':str(pair_lock),'reviewed':True,'results':[r]})
    return p


def test_replication_separates_scope_and_denominator(pair_lock,tmp_path):
    p=replication_manifest(pair_lock,tmp_path)
    result=evaluate(p,tmp_path/'w','rep');rs=list(read_tsv(result/'replication.tsv'))
    tested=next(r for r in rs if r['pair_id']=='insomnia__adhd')
    assert tested['status']=='REPLICATED' and tested['independence_scope']=='DISEASE_ONLY'
    assert float(tested['replication_family_fdr'])==pytest.approx(6e-10)
    assert len(rs)==6

@pytest.mark.parametrize('changes,expected',[
    ({'replication_cohorts':['D']},'NO_VALID_REPLICATION'),
    ({'disease_h2_z':2},'UNDERPOWERED'),
    ({'rg':-.4},'DISCORDANT'),
    ({'qc_pass':False},'NO_VALID_REPLICATION'),
    ({'independence_scope':'NONE'},'NO_VALID_REPLICATION'),
    ({'phenotype_match':'different'},'NO_VALID_REPLICATION')])
def test_replication_boundaries(pair_lock,tmp_path,changes,expected):
    p=replication_manifest(pair_lock,tmp_path,**changes)
    result=evaluate(p,tmp_path/'w','rep')
    assert next(r for r in read_tsv(result/'replication.tsv') if r['pair_id']=='insomnia__adhd')['status']==expected


def test_campaign_has_every_stage(pair_lock,tmp_path):
    p=tmp_path/'campaign.json';t=template(pair_lock,p)
    assert len(t['requirements'])==6*len(STAGES)
    result=audit(p,tmp_path/'w','audit');s=read_json(result/'status.json')
    assert s['requirements']==74 and s['resolved']==0 and not s['publication_ready']

@pytest.mark.parametrize('change',['remove','duplicate','fakepass'])
def test_campaign_cannot_hide_missing_work(pair_lock,tmp_path,change):
    p=tmp_path/'campaign.json';t=template(pair_lock,p);p.unlink()
    if change=='remove':t['requirements'].pop()
    elif change=='duplicate':t['requirements'].append(t['requirements'][0])
    else:t['requirements'][0]['status']='PASS'
    write_json(p,t)
    with pytest.raises(ContractError):audit(p,tmp_path/'w','audit')


def qtl_source(tmp_path,*,duplicate=False,complete=True):
    rows=[dict(zip(REQUIRED,['GENE1',1,100*i,'C','A',-.03,.01,.01,400,.7])) for i in range(1,11)]
    if duplicate:rows.append(rows[0])
    f=tmp_path/'qtl.tsv';write_tsv(f,REQUIRED,rows)
    source={'path':str(f),'sha256':sha256(f),'synthetic':True,'study_id':'SYNTHETIC',
            'genome_build':'GRCh37','ancestry':'EUR','molecular_trait':'expression',
            'tissue':'test cortex','source_uri':'synthetic://qtl','citation':'Artificial fixture',
            'n_semantics':'total','complete_cis_summary':complete,'column_map':{k:k for k in REQUIRED}}
    p=tmp_path/'qtl.json';write_json(p,source);return p


def test_qtl_duplicates_quarantined(tmp_path):
    p=qtl_source(tmp_path,duplicate=True);out=index(p,tmp_path/'w','qtl',synthetic=True)
    q=read_json(out/'qc.json');assert q['duplicated_coordinates']==1 and q['retained_rows']==9


def test_qtl_significant_only_rejected(tmp_path):
    p=qtl_source(tmp_path,complete=False)
    with pytest.raises(ContractError,match='Significant-only'):index(p,tmp_path/'w','qtl',synthetic=True)


def test_qtl_locus_full_pair_and_scope(source_factory,tmp_path):
    s,_,_=source_factory();g=normalize(s,tmp_path/'w','g',synthetic=True)
    q=index(qtl_source(tmp_path),tmp_path/'w','q',synthetic=True)
    p=locus(g,q,'GENE1',{'CHR':1,'START':1,'STOP':2000},tmp_path/'w','l',minimum_variants=2)
    rows=list(read_tsv(p/'pair.tsv.gz'))
    assert len(rows)==10 and all(float(r['BETA2'])==.03 for r in rows)
    assert all(float(r['N2'])==400 for r in rows)
    with pytest.raises(ContractError,match='Cis-only'):split_pair(p,tmp_path/'w','chunks')
    from brain6.ld import index_pair
    indexed=index_pair(p,tmp_path/'w','indexed')
    assert verify_artifact(indexed)['pair_metadata']['analysis_scope']=='GWAS_QTL_LOCUS'


def test_qtl_missing_feature_explicit(source_factory,tmp_path):
    s,_,_=source_factory();g=normalize(s,tmp_path/'w','g',synthetic=True)
    q=index(qtl_source(tmp_path),tmp_path/'w','q',synthetic=True)
    with pytest.raises(ContractError,match='COVERAGE'):
        locus(g,q,'ABSENT',{'CHR':1,'START':1,'STOP':2000},tmp_path/'w','l')

@pytest.mark.parametrize('bp,hit',[(99,False),(100,True),(101,True),(200,True),(201,False)])
def test_bed_coordinate_edges(tmp_path,bp,hit):
    p=tmp_path/'a.bed';p.write_text('chr1\t99\t200\nchr1\t150\t180\n')
    ix=IntervalIndex(p);assert ix.contains(1,bp)==hit
    assert len(ix.intervals[1])==1


def test_annotation_scores_are_pip_mass_not_double_counted(tmp_path):
    bed=tmp_path/'a.bed';bed.write_text('chr1\t99\t101\nchr1\t99\t101\n')
    a=tmp_path/'a.json';write_json(a,{'reviewed':True,'synthetic':True,'genome_build':'GRCh37',
        'independent_loci_reviewed':True,'matched_controls_reviewed':True,
        'annotations':[{'cell_type':'Cell1','path':str(bed),'sha256':sha256(bed),
        'genome_build':'GRCh37','coordinate_system':'BED0','citation':'Artificial','source_uri':'synthetic://bed'}]})
    rows=[dict(pair_id='insomnia__adhd',locus_id='L1',SNP='rs1',CHR=1,BP=100,PIP=.75,selected=1,stratum='s'),
          dict(pair_id='insomnia__adhd',locus_id='L1',SNP='rs2',CHR=1,BP=200,PIP=.25,selected=1,stratum='s')]
    v=tmp_path/'v.tsv';write_tsv(v,list(rows[0]),rows)
    out=score(v,a,tmp_path/'w','a',synthetic=True)
    assert float(next(read_tsv(out/'scores.tsv'))['score'])==.75


def test_new_python_dag_mr_resume(pair_lock,source_factory,tmp_path):
    from brain6.pipeline import run_pipeline
    s,_,_=source_factory()
    runtime=tmp_path/'runtime.json'
    write_json(runtime,{'schema_version':1,'pair_lock':str(pair_lock),'reviewed':True,'synthetic':True,
         'output_root':str(tmp_path/'dag'),'tasks':[
         {'id':'exposure','kind':'normalize','parameters':{'source':str(s)}},
         {'id':'mr_candidates','kind':'mr_candidates','depends_on':['exposure'],
          'parameters':{'exposure':'@exposure'}}]})
    first=run_pipeline(runtime,allow_synthetic=True)
    receipt=Path(first['mr_candidates'])/'receipt.json';digest=sha256(receipt)
    second=run_pipeline(runtime,allow_synthetic=True)
    assert first==second and sha256(receipt)==digest
    assert verify_artifact(first['mr_candidates'])['scientific_status']=='NO_SIGNAL'


def test_bed_sex_mito_and_alt_contigs_accounted(tmp_path):
    p=tmp_path/'a.bed';p.write_text('chrX\t0\t100\nchrM\t0\t100\nchr1_alt\t0\t100\nchr1\t0\t100\n')
    ix=IntervalIndex(p)
    assert ix.skipped_non_autosomal==3 and ix.contains(1,1)


@pytest.mark.parametrize('duplicate',[False,True])
def test_molecular_runner_blocks_without_native_statistics(pair_lock,source_factory,tmp_path,monkeypatch,duplicate):
    import numpy as np
    from brain6.molecular import run
    source,_,_=source_factory(trait_id='mdd')
    g=normalize(source,tmp_path/'data','gwas',synthetic=True)
    q=index(qtl_source(tmp_path),tmp_path/'data','qtl',synthetic=True)
    ld=tmp_path/'signed.tsv';np.savetxt(ld,np.eye(10),delimiter='\t')
    bim=tmp_path/'local.bim';bim.write_text(''.join(f'1 rs{i} 0 {100*i} A C\n' for i in range(1,11)))
    lm=tmp_path/'ld.json';write_json(lm,{'kind':'signed_r','counted_allele':'BIM_A1','synthetic':True,
        'ancestry':'EUR','genome_build':'GRCh37','matrix_path':str(ld),'matrix_sha256':sha256(ld),
        'bim_path':str(bim),'bim_sha256':sha256(bim)})
    query={'query_id':'mdd_GENE1_test','pair_id':'insomnia__mdd','gwas_trait':'mdd','feature_id':'GENE1',
        'gwas':str(g),'qtl_index':str(q),'region':{'CHR':1,'START':1,'STOP':2000},'ld_manifest':str(lm),
        'selection_reason':'Synthetic test, NOT a real gene','gwas_susie':{},'qtl_susie':{}}
    cfg={'reviewed':True,'feature_family_frozen':True,'pair_lock':str(pair_lock),
         'output_root':str(tmp_path/'molecular'),'queries':[query,query] if duplicate else [query],
         'minimum_variants':2,'max_working_bytes':2**27,'minimum_ld_coverage':.8,'susie_common':{},'coloc':{}}
    cp=tmp_path/'molecular.json';write_json(cp,cfg)
    def blocked(*args,**kwargs):
        raise ContractError('SYNTHETIC TEST: native engine intentionally unavailable')
    monkeypatch.setattr('brain6.molecular.run_job',blocked)
    if duplicate:
        with pytest.raises(ContractError,match='duplicate'):run(cp,allow_synthetic=True)
    else:
        r=run(cp,allow_synthetic=True)
        assert r['scientific_status']=='INSUFFICIENT_EVIDENCE'
        rows=list(read_tsv(tmp_path/'molecular/summary/molecular_evidence.tsv'))
        assert len(rows)==1 and rows[0]['status']=='FAILED_OR_BLOCKED' and rows[0]['max_h4']=='NA'
        assert verify_artifact(tmp_path/'molecular/mdd_GENE1_test_aligned')
        assert run(cp,allow_synthetic=True)['fingerprint']==r['fingerprint']


def test_mr_rejects_pleiotropy_clumping_policy(source_factory,tmp_path):
    c,l,y=mr_fixture(source_factory,tmp_path)
    with pytest.raises(ContractError,match='instrument independence'):
        harmonize(c,l,y,tmp_path/'w','bad',max_instrument_r2=.0001)
