import math
import copy
from pathlib import Path
import pytest
from brain6.io import *
from brain6.gwas import *
from brain6.artifacts import verify_artifact
from brain6.pairs import *

@pytest.mark.parametrize('a,b,expected',[
 (('A','C'),('A','C'),(1,'same')), (('A','C'),('C','A'),(-1,'swapped')),
 (('A','C'),('T','G'),(1,'complement')),(('A','C'),('G','T'),(-1,'complement_swapped')),
 (('A','T'),('A','T'),None),(('C','G'),('G','C'),None),(('A','C'),('A','G'),None),
 (('A','A'),('A','A'),None),(('AA','C'),('AA','C'),None)])
def test_alignment(a,b,expected):assert align(*a,*b)==expected

def test_normalize(source_factory,tmp_path):
    cfg,_,_=source_factory()
    path=normalize(cfg,tmp_path,'normal',synthetic=True)
    assert len(list(read_tsv(path/'sumstats.tsv.gz')))==10
    assert verify_artifact(path)['qc']['retained_rows']==10

def test_duplicate_all_removed(source_factory,tmp_path):
    _,_,rows=source_factory();rows.append(rows[0].copy())
    cfg,_,_=source_factory(rows)
    path=normalize(cfg,tmp_path,'dup',synthetic=True)
    rs=list(read_tsv(path/'sumstats.tsv.gz'))
    assert len(rs)==9 and 'rs1' not in {r['SNP'] for r in rs}

def test_ambiguous_positions_removed(source_factory,tmp_path):
    _,_,rows=source_factory();r=rows[0].copy();r['SNP']='alias1';rows.append(r)
    cfg,_,_=source_factory(rows);p=normalize(cfg,tmp_path,'ambiguous',synthetic=True)
    assert verify_artifact(p)['qc']['retained_rows']==9

@pytest.mark.parametrize('field,value,reason',[
 ('SE',0,'invalid_statistics'),('P',2,'invalid_statistics'),('N',1,'invalid_or_low_n'),
 ('A2','T','palindromic_or_non_biallelic_snp'),('EAF',.001,'frequency_filter'),
 ('INFO',.5,'info_filter'),('BP',-1,'non_autosomal_or_invalid_coordinate')])
def test_qc_filters(source_factory,tmp_path,field,value,reason):
    _,_,rows=source_factory();rows[0][field]=value
    cfg,_,_=source_factory(rows);p=normalize(cfg,tmp_path,'q',synthetic=True)
    qc=verify_artifact(p)['qc'];assert qc[reason]==1 and qc['retained_rows']==9

def test_or_conversion(source_factory,tmp_path):
    _,_,rows=source_factory();rows[0]['BETA']=2
    cfg,_,_=source_factory(rows,effect_scale='odds_ratio',se_scale='log_odds')
    p=normalize(cfg,tmp_path,'or',synthetic=True)
    r=next(read_tsv(p/'sumstats.tsv.gz'));assert float(r['BETA'])==pytest.approx(math.log(2))

def test_unreviewed_real_source_blocked(source_factory,tmp_path):
    cfg,_,_=source_factory(synthetic=False)
    with pytest.raises(ContractError):normalize(cfg,tmp_path,'no')

def test_checksum_blocks(source_factory,tmp_path):
    cfg,s,_=source_factory();Path(s['path']).write_text('corrupt')
    with pytest.raises(ContractError):normalize(cfg,tmp_path,'no',synthetic=True)

def test_source_malformed_fails_closed(source_factory,tmp_path):
    cfg,s,_=source_factory();p=Path(s['path']);p.write_text(p.read_text()+'bad\trow\n')
    s['sha256']=sha256(p);write_json(cfg,s,immutable=False)
    with pytest.raises(ContractError):normalize(cfg,tmp_path,'no',synthetic=True)
    assert not (tmp_path/'no').exists()

def test_join_swap(source_factory,tmp_path):
    cfg,_,rows=source_factory();left=normalize(cfg,tmp_path,'left',synthetic=True)
    for r in rows:r.update(A1='C',A2='A',BETA=-r['BETA'],EAF=.7)
    cfg,_,_=source_factory(rows);right=normalize(cfg,tmp_path,'right',synthetic=True)
    p=join_pair(left,right,tmp_path,'pair',min_variants=1,min_overlap=1,synthetic=True)
    rs=list(read_tsv(p/'pair.tsv.gz'));assert len(rs)==10
    assert all(float(r['BETA1'])==pytest.approx(float(r['BETA2'])) for r in rs)
    assert all(float(r['EAF2'])==pytest.approx(.3) for r in rs)

def test_build_mismatch(source_factory,tmp_path):
    cfg,_,_=source_factory();l=normalize(cfg,tmp_path,'l',synthetic=True)
    cfg,_,_=source_factory(genome_build='GRCh38');r=normalize(cfg,tmp_path,'r',synthetic=True)
    with pytest.raises(ContractError):join_pair(l,r,tmp_path,'pair',min_variants=1,synthetic=True)

def test_join_count_gate(source_factory,tmp_path):
    cfg,_,_=source_factory();l=normalize(cfg,tmp_path,'l',synthetic=True)
    cfg,_,_=source_factory();r=normalize(cfg,tmp_path,'r',synthetic=True)
    with pytest.raises(ContractError):join_pair(l,r,tmp_path,'no',min_variants=100,synthetic=True)
    assert not (tmp_path/'no').exists()

def test_rank_72(atlas_factory):
    m,c,_,_,cfg,_=atlas_factory();r=rank_pairs(m,cfg)
    assert len(r)==72 and sum(x['selection_rank']==1 for x in r)==6

def test_rank_does_not_change_fdr(atlas_factory):
    m,c,_,rows,cfg,_=atlas_factory();r=rank_pairs(m,cfg)
    assert r[0]['fdr']==rows[0]['fdr']

def test_missing_qc_no_guess(atlas_factory):
    m,c,_,_,cfg,_=atlas_factory();cfg['matrix_qc_column']='unknown'
    with pytest.raises(ContractError):rank_pairs(m,cfg)

def test_rank_duplicate_blocks(atlas_factory):
    m,c,_,rows,cfg,_=atlas_factory();rows[-1]=rows[0]
    write_tsv(m,list(rows[0]),rows,immutable=False)
    with pytest.raises(ContractError):rank_pairs(m,cfg)

def test_rank_invalid_not_selected(atlas_factory):
    m,c,_,rows,cfg,_=atlas_factory();rows[0]['rg']=1.3;rows[1]['analysis_status']='QC_FAILED_SENSITIVITY'
    write_tsv(m,list(rows[0]),rows,immutable=False);r=rank_pairs(m,cfg)
    assert not next(x for x in r if x['pair_id']=='insomnia__adhd')['eligible']
    assert not next(x for x in r if x['pair_id']=='sleepdur__adhd')['eligible']

def test_freeze_threshold_and_tampering(atlas_factory,tmp_path):
    m,c,d,_,_,_=atlas_factory();out=tmp_path/'lock.json'
    lock=freeze_pairs(m,c,d,out,'test-reviewer',synthetic=True)
    assert lock['genomewide_family_threshold']==pytest.approx(5e-8/6)
    assert verify_lock(out)['n_pairs']==6
    lock['n_pairs']=2;write_json(out,lock,immutable=False)
    with pytest.raises(ContractError):verify_lock(out)

def test_unreviewed_pair_not_freezable(atlas_factory,tmp_path):
    m,c,d,rows,cfg,decisions=atlas_factory();decisions[1]['sleep_trait']='imaginary'
    write_tsv(d,list(decisions[0]),decisions,immutable=False)
    with pytest.raises(ContractError):freeze_pairs(m,c,d,tmp_path/'lock','test',synthetic=True)

def test_synthetic_lock_cannot_be_empirical(atlas_factory,tmp_path):
    m,c,d,*_=atlas_factory()
    with pytest.raises(ContractError):freeze_pairs(m,c,d,tmp_path/'l','review')
