"""Execute actual pytest suites and retain JUnit/log evidence without hiding skips."""
from __future__ import annotations
import os
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path
from .io import require,read_json,write_json,write_tsv
from .artifacts import transaction
from .doctor import inspect_environment


def junit_counts(path):
    tree=ET.parse(path)
    cases=tree.findall('.//testcase');bad=tree.findall('.//error')
    require(cases or bad,'Empty JUnit output is not a passing test run')
    counts={'tests':len(cases),'failures':len(tree.findall('.//failure')),
            'errors':len(bad),'skipped':len(tree.findall('.//skipped'))}
    counts['passed']=len(cases)-counts['failures']-counts['errors']-counts['skipped']
    return counts


def run(manifest,root,name):
    c=read_json(manifest);require(c.get('reviewed') is True,'Freeze test suite scope before execution')
    mode=c.get('mode','regression');require(mode in {'regression','native_environment'},'Unknown validation mode')
    suites=c['test_suites'];require(suites and len(suites)==len(set(s['id'] for s in suites)),'Empty/duplicate test suite scope')
    inputs=[manifest];paths=[]
    for s in suites:
        p=Path(s['path']).resolve();require(p.exists(),'Missing declared test suite')
        files=sorted(p.rglob('*.py')) if p.is_dir() else [p]
        require(files,'No test source files');inputs.extend(files);paths.append(p)
    if mode=='native_environment':
        canonical=Path(__file__).resolve().parents[1]/'tests/test_native_runtime.py'
        require(paths==[canonical.resolve()],'Native validation must execute the shipped actual-runtime tests')
    stage='native_environment_validation' if mode=='native_environment' else 'regression_validation'
    timeout=int(c.get('timeout_seconds',1800));require(timeout>=1,'Invalid validation timeout')
    with transaction(root,name,stage=stage,inputs=inputs,parameters=c,synthetic=c.get('synthetic',False)) as (work,meta):
        environment=inspect_environment(work/'environment.json')
        outcomes=[]
        for suite,p in zip(suites,paths):
            from .io import safe_id
            sid=safe_id(suite['id']);xml=work/(sid+'.junit.xml');log=work/(sid+'.log')
            argv=[sys.executable,'-m','pytest','-q','--disable-warnings',str(p),'--junitxml='+str(xml)]
            with log.open('w') as f:
                try:
                    result=subprocess.run(argv,stdout=f,stderr=subprocess.STDOUT,timeout=timeout,check=False,
                       cwd=str(c.get('working_directory',p.parent)),env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1'})
                    code=result.returncode
                except subprocess.TimeoutExpired:code=124
            counts=junit_counts(xml) if xml.is_file() else {'tests':0,'passed':0,'failures':0,'errors':1,'skipped':0}
            outcomes.append({'suite':sid,'exit_code':code,**counts})
        good=all(r['exit_code']==0 and r['tests']>0 and r['failures']==r['errors']==r['skipped']==0 for r in outcomes)
        if mode=='native_environment':good=good and not environment['missing_requirements']
        meta['scientific_status']='PASS' if good else 'INSUFFICIENT_EVIDENCE'
        write_tsv(work/'test_suites.tsv',list(outcomes[0]),outcomes)
        write_json(work/'status.json',{'status':meta['scientific_status'],'mode':mode,
           'tests':sum(r['tests'] for r in outcomes),'passed':sum(r['passed'] for r in outcomes),
           'skipped':sum(r['skipped'] for r in outcomes),'native_functional_scope':['R syntax','SuSiE/coloc artificial data','PLINK artificial LD/clumping','TwoSampleMR artificial data'] if mode=='native_environment' else [],
           'note':'Tests are software validation only, not real-GWAS results. Skips never become executed analyses. Native functional tests do not prove all production models valid.'})
    return Path(root)/name
