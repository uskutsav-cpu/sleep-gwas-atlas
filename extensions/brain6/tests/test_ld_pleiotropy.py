import json
from pathlib import Path
import numpy as np
import pytest
from brain6.io import *
from brain6.gwas import normalize,join_pair
from brain6.ld import *
from brain6.artifacts import verify_artifact
from brain6.pleiotropy import split_pair,collate_placo
from brain6.stats import bh

@pytest.mark.parametrize('matrix',[
 [[1,.4],[.4,1]],[[1,-.5],[-.5,1]],[[1,1],[1,1]]])
def test_ld_valid(matrix):assert validate_ld(matrix)>=-1e-6

@pytest.mark.parametrize('matrix',[
 [[1,2],[2,1]],[[1,.2],[.3,1]],[[1,float('nan')],[float('nan'),1]],
 [[.9,0],[0,1]],[[1,1,1],[1,1,-1],[1,-1,1]],[[1]]])
def test_ld_invalid(matrix):
    with pytest.raises(ContractError):validate_ld(matrix)

def test_block_mapping_no_double_count(tmp_path):
    b=tmp_path/'blocks.tsv';write_tsv(b,['LOC','CHR','START','STOP'],[
        {'LOC':'L1','CHR':1,'START':1,'STOP':1000},{'LOC':'L2','CHR':1,'START':1001,'STOP':2000}])
    l=tmp_path/'leads.tsv';write_tsv(l,['SNP','CHR','BP','P'],[
        {'SNP':'rs1','CHR':1,'BP':100,'P':1e-10},{'SNP':'rs2','CHR':1,'BP':900,'P':1e-12}])
    rows=map_loci(l,b,tmp_path/'loci.tsv')
    assert len(rows)==1 and rows[0]['n_clump_leads']==2 and rows[0]['lead_snp']=='rs2'

def test_overlap_blocks_rejected(tmp_path):
    p=tmp_path/'b';write_tsv(p,['LOC','CHR','START','STOP'],[
        {'LOC':'a','CHR':1,'START':1,'STOP':100},{'LOC':'b','CHR':1,'START':100,'STOP':200}])
    with pytest.raises(ContractError):read_blocks(p)

def test_ld_aligned_both_traits(source_factory,tmp_path):
    cfg,_,_=source_factory();l=normalize(cfg,tmp_path,'l',synthetic=True)
    cfg,_,_=source_factory();r=normalize(cfg,tmp_path,'r',synthetic=True)
    pair=join_pair(l,r,tmp_path,'p',min_variants=1,synthetic=True)
    idx=index_pair(pair,tmp_path,'idx')
    matrix=tmp_path/'ld.tsv';np.savetxt(matrix,np.eye(10),delimiter='\t')
    bim=tmp_path/'locus.bim'
    bim.write_text(''.join(f'1 rs{j} 0 {j*100} C A\n' for j in range(1,11)))
    manifest={'kind':'signed_r','counted_allele':'BIM_A1','matrix_path':str(matrix),'matrix_sha256':sha256(matrix),
        'bim_path':str(bim),'bim_sha256':sha256(bim),'genome_build':'GRCh37','ancestry':'EUR','synthetic':True}
    mp=tmp_path/'ld.json';write_json(mp,manifest)
    out=prepare_locus(idx,mp,{'CHR':1,'START':1,'STOP':1000,'lead_snp':'rs1'},tmp_path,'aligned',synthetic=True)
    rows=list(read_tsv(out/'locus.tsv'))
    assert all(float(x['BETA1'])==pytest.approx(-.1) and float(x['BETA2'])==pytest.approx(-.1) for x in rows)
    assert all(float(x['EAF1'])==pytest.approx(.7) for x in rows)
    manifest['kind']='r_squared';write_json(mp,manifest,immutable=False)
    with pytest.raises(ContractError):prepare_locus(idx,mp,{'CHR':1,'START':1,'STOP':1000},tmp_path,'bad',synthetic=True)

def test_split_preserves_rows(source_factory,tmp_path):
    cfg,_,_=source_factory();l=normalize(cfg,tmp_path,'l',synthetic=True)
    cfg,_,_=source_factory();r=normalize(cfg,tmp_path,'r',synthetic=True)
    pair=join_pair(l,r,tmp_path,'p',min_variants=1,synthetic=True)
    split=split_pair(pair,tmp_path,'parts',chunk_rows=3)
    m=read_json(split/'chunks.json');assert [r['rows'] for r in m['chunks']]==[3,3,3,1]

def placo_fixture(tmp_path,fail=False):
    fields=['SNP','CHR','BP','Z1','Z2','P_PLACO','status']
    rows=[]
    ps=[.01,.01,.04,.8,1e-12,None]
    for i,p in enumerate(ps):
        rows.append(dict(zip(fields,[f'rs{i}',1,100+i,1,1,'NA' if p is None else p,
                   'NUMERICAL_FAILURE' if p is None and fail else 'EXCLUDED_EXTREME_Z' if p is None else 'TESTED'])))
    files=[]
    for i,chunk in enumerate([rows[:3],rows[3:]]):
        path=tmp_path/f'chunk{i}.tsv';write_tsv(path,fields,chunk)
        files.append({'path':str(path),'sha256':sha256(path)})
    m=tmp_path/'manifest.tsv';write_tsv(m,['path','sha256'],files)
    return m,ps

def test_global_bh_not_chunk_bh(tmp_path):
    manifest,ps=placo_fixture(tmp_path)
    out=collate_placo(manifest,tmp_path,'global',expected_rows=6,n_pairs=6,synthetic=True)
    rows=list(read_tsv(out/'results.tsv.gz'))
    expected=bh(ps[:5])
    assert np.allclose([float(r['Q_WITHIN_PAIR']) for r in rows[:5]],expected)
    status=read_json(out/'status.json')
    assert status['denominator']==5 and status['headline_variants']==1
    assert rows[-1]['Q_WITHIN_PAIR']=='NA'

def test_failures_remain_in_denominator(tmp_path):
    manifest,ps=placo_fixture(tmp_path,fail=True)
    out=collate_placo(manifest,tmp_path,'global',expected_rows=6,n_pairs=6,max_failure_rate=.01,synthetic=True)
    rows=list(read_tsv(out/'results.tsv.gz'))
    assert np.allclose([float(r['Q_WITHIN_PAIR']) for r in rows[:5]],bh(ps)[:5])
    assert verify_artifact(out)['scientific_status']=='FAILED_QC_NOT_CONSUMED'
    assert rows[-1]['status']=='NUMERICAL_FAILURE' and rows[-1]['P_PLACO']=='NA'

def test_missing_chunk_cannot_finish(tmp_path):
    manifest,_=placo_fixture(tmp_path)
    with pytest.raises(ContractError):collate_placo(manifest,tmp_path,'global',expected_rows=7,n_pairs=6,synthetic=True)
    assert not (tmp_path/'global').exists()
