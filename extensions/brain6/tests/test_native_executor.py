import io
import json
import sys
from pathlib import Path
from unittest.mock import patch
import pytest
from brain6.io import *
from brain6.native import parse_h2,parse_rg
from brain6.executor import validate_job,run_job
from brain6.pipeline import sorted_tasks,expand
from brain6.download import download


def test_parse_h2():
    r=parse_h2('Total Observed scale h2: 0.123 (0.02)\nIntercept: 1.04 (0.01)\n')
    assert r['z']==pytest.approx(6.15) and r['scale']=='Observed'

def test_parse_liability_h2_not_capped():
    # Parser preserves the value; adapter's QC rejects >1 instead of hiding it.
    r=parse_h2('Total Liability scale h2: 5.98 (0.1)\nIntercept: 1.1 (0.01)\n')
    assert r['h2']==5.98

def test_parse_h2_failure():
    with pytest.raises(ContractError):parse_h2('some partial log')

def test_parse_rg():
    text='Summary of Genetic Correlation Results\np1 p2 rg se z p gcov_int\na.gz b.gz 0.3 0.05 6 1e-9 0.02\n\nAnalysis finished'
    r=parse_rg(text);assert len(r)==1 and r[0]['rg']==.3

@pytest.mark.parametrize('rg',['nan','1.4','-1.1'])
def test_rg_nonphysical(rg):
    with pytest.raises(ContractError):parse_rg(f'p1 p2 rg se p\na b {rg} 0.1 0.01\n')


def fixture_job(tmp_path,status='PASS'):
    script=tmp_path/'native_stub.py'
    script.write_text("import os,json\nfrom pathlib import Path\np=Path(os.environ['BRAIN6_JOB_DIR'])\n(p/'result.tsv').write_text('test\\n')\n(p/'status.json').write_text(json.dumps({'status':"+repr(status)+"}))\n")
    job={'schema_version':1,'job_id':'stub','reviewed':True,'synthetic':True,
         'argv':[sys.executable,'{input:script}'],'inputs':{'script':file_record(script)},
         'outputs':{'result':'result.tsv','semantic_status':'status.json'},'minimum_free_bytes':0,'threads':1}
    path=tmp_path/'job.json';write_json(path,job)
    return path,job


def test_executor_runs_and_resumes(tmp_path):
    path,_=fixture_job(tmp_path);a=run_job(path,tmp_path/'run',allow_synthetic=True)
    b=run_job(path,tmp_path/'run',allow_synthetic=True)
    assert a['fingerprint']==b['fingerprint'] and a['scientific_status']=='PASS'


def test_executor_synthetic_guard(tmp_path):
    path,_=fixture_job(tmp_path)
    with pytest.raises(ContractError):run_job(path,tmp_path/'run')


def test_executor_changed_input_refused(tmp_path):
    path,job=fixture_job(tmp_path);Path(job['inputs']['script']['path']).write_text('changed')
    with pytest.raises(ContractError):run_job(path,tmp_path/'run',allow_synthetic=True)


def test_executor_semantic_failure(tmp_path):
    path,_=fixture_job(tmp_path,status='FAILED_QC')
    with pytest.raises(ContractError):run_job(path,tmp_path/'run',allow_synthetic=True)
    assert not (tmp_path/'run/stub').exists()


def test_executor_low_power_not_null(tmp_path):
    path,_=fixture_job(tmp_path,status='INSUFFICIENT_EVIDENCE')
    a=run_job(path,tmp_path/'run',allow_synthetic=True)
    assert a['scientific_status']=='INSUFFICIENT_EVIDENCE'


def test_shell_jobs_disallowed(tmp_path):
    _,j=fixture_job(tmp_path);j['argv']=['bash','-c','echo unsafe']
    with pytest.raises(ContractError):validate_job(j)


def test_dag_sort_and_cycle():
    tasks=[{'id':'b','depends_on':['a']},{'id':'a','depends_on':[]}]
    assert [t['id'] for t in sorted_tasks(tasks)]==['a','b']
    tasks[1]['depends_on']=['b']
    with pytest.raises(ContractError):sorted_tasks(tasks)


def test_expansion_safe(tmp_path):
    assert expand('@job/a.tsv',{'job':tmp_path})==str(tmp_path/'a.tsv')
    with pytest.raises(ContractError):expand('@job/../../x',{'job':tmp_path})

class FakeResponse(io.BytesIO):
    def __init__(self,body,status=200,headers=None):
        super().__init__(body);self.status=status;self.headers=headers or {}
    def geturl(self):return 'https://example.org/data'


def test_download_verified(tmp_path):
    import hashlib
    body=b'test-data';h=hashlib.sha256(body).hexdigest()
    with patch('urllib.request.urlopen',return_value=FakeResponse(body,headers={'Content-Length':str(len(body))})):
        path=download('https://example.org/data',tmp_path/'download',h,len(body))
    assert path.read_bytes()==body


def test_resume_range_is_checked(tmp_path):
    import hashlib
    body=b'1234567890';h=hashlib.sha256(body).hexdigest();d=tmp_path/'d'
    (tmp_path/'d.partial').write_bytes(body[:4])
    write_json(tmp_path/'d.download.json',{'url':'https://example.org/data','sha256':h,'bytes':10})
    with patch('urllib.request.urlopen',return_value=FakeResponse(body[4:],206,{'Content-Range':'bytes 4-9/10','Content-Length':'6'})):
        download('https://example.org/data',d,h,10)
    assert d.read_bytes()==body


def test_ignored_range_preserves_partial(tmp_path):
    import hashlib
    body=b'1234567890';h=hashlib.sha256(body).hexdigest();d=tmp_path/'d'
    (tmp_path/'d.partial').write_bytes(body[:4])
    write_json(tmp_path/'d.download.json',{'url':'https://example.org/data','sha256':h,'bytes':10})
    with patch('urllib.request.urlopen',return_value=FakeResponse(body,200)):
        with pytest.raises(ContractError):download('https://example.org/data',d,h,10)
    assert (tmp_path/'d.partial').read_bytes()==body[:4] and not d.exists()


def test_download_bad_hash_not_promoted(tmp_path):
    with patch('urllib.request.urlopen',return_value=FakeResponse(b'123')):
        with pytest.raises(ContractError):download('https://example.org/data',tmp_path/'d','0'*64,3)
    assert not (tmp_path/'d').exists()


def test_download_http_forbidden(tmp_path):
    with pytest.raises(ContractError):download('http://example.org/data',tmp_path/'d','0'*64,3)
