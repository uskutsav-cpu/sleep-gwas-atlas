"""Read-only environment audit. Never install, guess versions, or equate import with validation."""
from __future__ import annotations
import importlib.metadata
import json
import os
import platform
import shutil
import subprocess
import sys
from pathlib import Path
from .io import file_record, json_hash, write_json, require

R_PACKAGES=('R.utils','data.table','jsonlite','susieR','coloc','LAVA','GenomicSEM','TwoSampleMR')

def inspect_environment(out=None, *, expected=None):
    result={'schema_version':1,'python':{'executable':sys.executable,'version':platform.python_version()},
            'platform':platform.platform(),'machine':platform.machine(),'tools':{},'python_packages':{},
            'r_packages':{},'native_statistics_executed':False}
    for p in ('numpy','scipy','matplotlib','pytest','snakemake'):
        try: result['python_packages'][p]=importlib.metadata.version(p)
        except importlib.metadata.PackageNotFoundError:result['python_packages'][p]=None
    for cmd,args in [('git',['--version']),('Rscript',['--version']),('plink',['--version']),
                     ('plink1.9',['--version']),('plink2',['--version']),('snakemake',['--version']),
                     ('matlab',None)]:
        path=shutil.which(cmd)
        row={'path':path,'status':'NOT_FOUND' if path is None else 'FOUND_NOT_STATISTICALLY_VALIDATED'}
        if path:
            row['binary_sha256']=file_record(Path(path).resolve())['sha256']
            if args:
                try:
                    p=subprocess.run([path,*args],capture_output=True,text=True,timeout=20,check=False)
                    row.update(returncode=p.returncode,version_output=(p.stdout+p.stderr).strip()[:1500])
                except (OSError,subprocess.TimeoutExpired) as e:row.update(status='PROBE_FAILED',error=str(e))
        result['tools'][cmd]=row
    r=result['tools']['Rscript']['path']
    for p in R_PACKAGES:result['r_packages'][p]={'status':'NOT_PROBED_R_ABSENT'}
    if r:
        # Base-R TSV output works even if jsonlite itself is absent.
        code=';'.join(["p<-c("+','.join(json.dumps(x) for x in R_PACKAGES)+")",
             "for(x in p){ok<-requireNamespace(x,quietly=TRUE);cat(x,if(ok)as.character(packageVersion(x))else'ABSENT',sep='\\t');cat('\\n')}"])
        try:
            p=subprocess.run([r,'--vanilla','-e',code],capture_output=True,text=True,timeout=90,check=False)
            require(p.returncode==0,'R package probe failed: '+p.stderr[:1000])
            found={line.split('\t')[0]:line.split('\t')[1] for line in p.stdout.splitlines() if '\t' in line}
            require(set(found)==set(R_PACKAGES),'Incomplete native package probe')
            result['r_packages']={k:{'status':'ABSENT' if v=='ABSENT' else 'IMPORTABLE_NOT_STATISTICALLY_VALIDATED',
                                      'version':None if v=='ABSENT' else v} for k,v in found.items()}
        except (OSError,subprocess.TimeoutExpired,ValueError) as e:result['r_probe_error']=str(e)
    result['missing_requirements']=[k for k,v in result['r_packages'].items() if not v.get('version')]
    if not any(result['tools'][k]['path'] for k in ('plink','plink1.9')):
        result['missing_requirements'].append('PLINK_1.9')
    result['fingerprint']=json_hash({k:v for k,v in result.items() if k!='fingerprint'})
    if expected:
        from .io import read_json
        require(read_json(expected)['fingerprint']==result['fingerprint'],'Runtime changed from pinned environment')
    if out:write_json(out,result)
    return result
