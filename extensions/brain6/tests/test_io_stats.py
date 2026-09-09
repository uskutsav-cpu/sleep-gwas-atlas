import gzip
import math
from pathlib import Path
import numpy as np
import pytest
from scipy.stats import false_discovery_control
from brain6.io import *
from brain6.artifacts import transaction,verify_artifact
from brain6.stats import *

@pytest.mark.parametrize('bad',['../x','a/../../b','/etc/passwd'])
def test_traversal(tmp_path,bad):
    with pytest.raises(ContractError):safe_write_path(tmp_path,bad)

@pytest.mark.parametrize('bad',['../x','a/b','','x y',';rm','..'])
def test_bad_id(bad):
    with pytest.raises(ContractError):safe_id(bad)

def test_symlink_escape(tmp_path):
    (tmp_path/'link').symlink_to('/tmp')
    with pytest.raises(ContractError):safe_write_path(tmp_path,'link/result')

def test_atomic_abort(tmp_path):
    p=tmp_path/'result.tsv'
    with pytest.raises(RuntimeError):
        with atomic_text(p) as f:
            f.write('half');raise RuntimeError()
    assert not p.exists()
    assert not list(tmp_path.glob('*.partial'))

def test_atomic_no_clobber(tmp_path):
    write_json(tmp_path/'a.json',{'x':1})
    with pytest.raises(ContractError):write_json(tmp_path/'a.json',{'x':2})
    assert read_json(tmp_path/'a.json')=={'x':1}

def test_deterministic_gzip_and_magic(tmp_path):
    for name in ['a.gz','b.gz']:
        with atomic_text(tmp_path/name) as f:f.write('SNP\tP\nrs1\t0.1\n')
    assert sha256(tmp_path/'a.gz')==sha256(tmp_path/'b.gz')
    (tmp_path/'b.gz').rename(tmp_path/'not_gzip.txt')
    assert list(read_tsv(tmp_path/'not_gzip.txt'))[0]['SNP']=='rs1'

def test_malformed_tsv(tmp_path):
    p=tmp_path/'t';p.write_text('a\tb\n1\t2\t3\n')
    with pytest.raises(ContractError):list(read_tsv(p))

def test_transaction_corruption(tmp_path):
    inp=tmp_path/'input';inp.write_text('verified')
    with transaction(tmp_path,'stage',stage='test',inputs=[inp],parameters={},synthetic=True) as (w,m):
        (w/'x').write_text('complete')
    assert verify_artifact(tmp_path/'stage')['synthetic']
    (tmp_path/'stage/x').write_text('corrupt')
    with pytest.raises(ContractError):verify_artifact(tmp_path/'stage')

def test_lock_exclusive(tmp_path):
    p=tmp_path/'lock'
    with exclusive_lock(p):
        with pytest.raises(ContractError):
            with exclusive_lock(p):pass
    assert not p.exists()

@pytest.mark.parametrize('p',[[.04,.01,.03,.001],[.01,.01,.05,1],[0,1,1],[.02]])
def test_bh_matches_scipy(p):
    assert np.allclose(bh(p),false_discovery_control(np.asarray(p,dtype=float)))

def test_bh_missing_denominator():
    q=bh([.01,None,.04,float('nan')])
    assert q==[.04,None,.08,None]
    assert bh([.01,.04],family_size=4)==[.04,.08]

@pytest.mark.parametrize('bad',[-1,1.1,float('inf')])
def test_bh_invalid(bad):
    with pytest.raises(ContractError):bh([bad])

def test_bh_bad_family():
    with pytest.raises(ContractError):bh([.1,.2],1)

def test_effective_n():assert effective_n(100,100)==200

def test_ivw_known_effect():
    x=np.array([.1,.2,.15,.3]);ans=ivw(x,[.01]*4,2*x,[.02]*4)
    assert ans['beta']==pytest.approx(2)
    assert ans['Q']==pytest.approx(0,abs=1e-12)
    assert ans['min_F']==pytest.approx(100)

@pytest.mark.parametrize('ind,match,p,rg,hz,expected',[
 (True,'exact',.0001,.3,5,'REPLICATED'),(False,'exact',.0001,.3,5,'NO_VALID_REPLICATION'),
 (True,'exact',.1,.3,5,'DIRECTIONALLY_SUPPORTED'),(True,'exact',.001,-.3,5,'DISCORDANT'),
 (True,'exact',.1,-.3,5,'INCONCLUSIVE_OPPOSITE_SIGN'),(True,'exact',.001,.3,2,'UNDERPOWERED')])
def test_replication(ind,match,p,rg,hz,expected):
    assert compare_replication(.3,rg,.1,p,alpha=.01,h2_z=hz,independent=ind,phenotype_match=match)==expected

def test_permutation_reproducible():
    score=np.arange(20);selected=np.array([False]*10+[True]*10);strata=['x']*20
    a=stratified_enrichment(score,selected,strata,permutations=199,seed=1)
    b=stratified_enrichment(score,selected,strata,permutations=199,seed=1)
    assert a==b and 0<a['p']<=.02

def test_permutation_no_controls():
    with pytest.raises(ContractError):stratified_enrichment([1,2],[True,False],['x','y'],permutations=99)


@pytest.mark.parametrize('field,value',[('discovery_rg',float('nan')),('discovery_rg',2),('h2_z',float('nan')),('alpha',0),('alpha',2)])
def test_replication_rejects_invalid_policy(field,value):
    from brain6.stats import compare_replication
    args=dict(discovery_rg=.3,replication_rg=.3,replication_se=.05,replication_p=.001,
              alpha=.01,h2_z=5,independent=True,phenotype_match="exact")
    args[field]=value
    with pytest.raises(ContractError):compare_replication(**args)
