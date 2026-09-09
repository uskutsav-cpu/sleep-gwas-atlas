"""Actual native programs on artificial data, never mocked.

Absent runtimes/packages SKIP explicitly. Passing Python contracts alone do not
imply these have run. All input data in this file are synthetic, not GWAS data.
"""
import json
import os
from pathlib import Path
import shutil
import subprocess
import numpy as np
import pytest
from brain6.io import write_json, write_tsv, read_json, read_tsv
from brain6.pipeline import create_job
from brain6.executor import run_job

pytestmark=pytest.mark.native
ROOT=Path(__file__).resolve().parents[1]


def native_r(packages=()):
    r=shutil.which('Rscript')
    if not r:pytest.skip('Rscript absent: native statistics NOT EXECUTED')
    for package in packages:
        probe=subprocess.run([r,'-e',f'quit(status=if(requireNamespace("{package}",quietly=TRUE)) 0 else 2)'],capture_output=True,text=True)
        if probe.returncode:pytest.skip(f'Native R package absent: {package}; NOT EXECUTED')
    return r


def adapter(tmp_path,method,cfg,job_name):
    settings=tmp_path/f'{job_name}.json';write_json(settings,cfg)
    job=create_job(settings,method=method,job_id=job_name,inputs={},
        outputs={'semantic_status':'status.json'},root=tmp_path/'jobs',
        reviewed=True,synthetic=True,minimum_free_bytes=0)
    run_job(job,tmp_path/'native',allow_synthetic=True)
    return tmp_path/'native'/job_name


def test_real_r_parse_all_adapters():
    r=native_r()
    code='for (f in commandArgs(TRUE)) { parse(file=f); cat("PARSED",basename(f),"\\n") }'
    result=subprocess.run([r,'-e',code,*map(str,sorted(set((ROOT/'brain6/resources').glob('*.R'))|set((ROOT/'scripts').glob('*.R'))))],capture_output=True,text=True)
    assert result.returncode==0,result.stdout+result.stderr


def test_real_r_susie_coloc_synthetic(tmp_path):
    native_r(('data.table','jsonlite','susieR','coloc'))
    rng=np.random.default_rng(92704);n=10000;p=40
    X=rng.binomial(2,.3,size=(1200,p)).astype(float)
    R=np.corrcoef(X,rowvar=False)
    z1=14*R[:,19]+rng.multivariate_normal(np.zeros(p),R)
    z2=13*R[:,19]+rng.multivariate_normal(np.zeros(p),R)
    ld=tmp_path/'signed.tsv';np.savetxt(ld,R,delimiter='\t',fmt='%.15g')
    locus=tmp_path/'locus.tsv'
    rows=[dict(SNP=f'rs{i+1}',CHR=1,BP=100000+i*100,A1='A',A2='C',
        BETA1=z1[i]/np.sqrt(n),SE1=1/np.sqrt(n),N1=n,
        BETA2=z2[i]/np.sqrt(n),SE2=1/np.sqrt(n),N2=n) for i in range(p)]
    write_tsv(locus,list(rows[0]),rows)
    cfg=dict(locus_file=str(locus),ld_file=str(ld),trait_type='quant',sample_size=n,
        sample_size_justification='Known synthetic N, not a real study',sdY=1,
        max_relative_n_deviation=.01,max_ld_mismatch=.3,L=3,coverage=.95,min_abs_corr=.5,
        max_iterations=500,seed=92704)
    a=adapter(tmp_path,'susie',{**cfg,'trait_number':1},'fit1')
    b=adapter(tmp_path,'susie',{**cfg,'trait_number':2},'fit2')
    a,b=Path(a),Path(b)
    assert all(0<=float(r['PIP'])<=1 for r in read_tsv(a/'pip.tsv'))
    out=Path(adapter(tmp_path,'coloc',dict(fit1=str(a/'fit.rds'),fit2=str(b/'fit.rds'),
        p1=1e-4,p2=1e-4,p12=1e-5,p12_sensitivity=[1e-6,1e-5],seed=92704),'coloc'))
    assert read_json(out/'status.json')['status'] in {'PASS','NO_SIGNAL','INSUFFICIENT_EVIDENCE'}
    # Model outputs, not pre-imposed high H4 or a real shared-signal claim.


def test_real_plink_clump_and_signed_ld(tmp_path):
    plink=shutil.which('plink') or shutil.which('plink1.9')
    if not plink:pytest.skip('PLINK 1.9 absent: native LD/clumping NOT EXECUTED')
    version=subprocess.run([plink,'--version'],capture_output=True,text=True)
    if '1.9' not in version.stdout+version.stderr:pytest.skip('This adapter requires PLINK 1.9')
    rng=np.random.default_rng(582);p=12
    prefix=tmp_path/'tiny'
    prefix.with_suffix('.map').write_text(''.join(f'1 rs{i+1} 0 {10000+i*100}\n' for i in range(p)))
    with prefix.with_suffix('.ped').open('w') as f:
        for i in range(100):
            g=rng.binomial(2,.3,size=p);alleles=[]
            for count in g:alleles.extend([('A','A'),('A','C'),('C','C')][count])
            f.write(f'F{i} I{i} 0 0 1 -9 '+' '.join(alleles)+'\n')
    result=subprocess.run([plink,'--file',str(prefix),'--make-bed','--out',str(prefix)],capture_output=True,text=True)
    assert result.returncode==0,result.stderr
    assoc=tmp_path/'assoc.tsv';write_tsv(assoc,['SNP','P'],[dict(SNP=f'rs{i+1}',P=1e-9 if i<3 else .5) for i in range(p)])
    out=Path(adapter(tmp_path,'clump',dict(mode='clump',plink=plink,reference_prefix=str(prefix),
        associations=str(assoc),p_column='P',p1=5e-8,p2=1,r2=.1,kb=500,threads=1),'clump'))
    assert len(list(read_tsv(out/'leads.tsv')))>=1
    ld=Path(adapter(tmp_path,'ld',dict(mode='ld',plink=plink,reference_prefix=str(prefix),
        reference_ancestry='EUR',genome_build='GRCh37',chromosome=1,start=10000,stop=12000,
        max_working_bytes=2**28,threads=1),'ld'))
    r=np.loadtxt(ld/'signed.ld')
    assert r.shape==(p,p) and np.allclose(r,r.T,atol=1e-6)


def test_real_twosamplemr_synthetic(tmp_path):
    native_r(('data.table','jsonlite','TwoSampleMR'))
    rng=np.random.default_rng(182);x=[];y=[]
    for i in range(18):
        b=.2+.012*i
        row=dict(SNP=f'rs{i+1}',BETA=b,SE=.02,P=1e-14,A1='A',A2='C',EAF=.3,N=100000)
        x.append(row);y.append({**row,'BETA':.4*b+rng.normal(0,.018),'P':.1})
    ex=tmp_path/'ex.tsv';oy=tmp_path/'oy.tsv';write_tsv(ex,list(x[0]),x);write_tsv(oy,list(y[0]),y)
    out=Path(adapter(tmp_path,'mr',dict(exposure=str(ex),outcome=str(oy),instrument_p=5e-8,minimum_F=10,
        instruments_ld_clumped=True,instrument_selection_exposure_only=True,sample_overlap_reviewed=True,seed=182),'mr'))
    result=list(read_tsv(out/'mr.tsv'))
    assert sum(r['method']=='Inverse variance weighted' for r in result)==1
