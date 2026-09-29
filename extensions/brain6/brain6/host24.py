"""One-command read-only host/snapshot investigation; no production data download or jobs."""
from __future__ import annotations
import subprocess
from pathlib import Path
from .io import require, write_json, write_tsv, read_json, file_record
from .doctor import inspect_environment
from .atlas_review import review, inspect_covariance, review_replication
from .inputs24 import inventory
from .completion24 import STEPS


def inspect(repo, out, *, archive_date=None, hash_dense=False, scan_shards=False):
    repo=Path(repo).resolve();out=Path(out).resolve()
    require(repo.is_dir(),'Repository/snapshot does not exist')
    require(not out.exists(),'Choose a new immutable audit output directory')
    require(not out.is_relative_to(repo/'data') and not out.is_relative_to(repo/'ref'),
            'Do not write audit outputs into raw/reference data')
    out.mkdir(parents=True)
    env=inspect_environment(out/'environment.json')
    git=subprocess.run(['git','-C',str(repo),'rev-parse','HEAD'],capture_output=True,text=True,check=False)
    state={'schema_version':1,'source_root':str(repo),'archive_date':archive_date,
       'git_commit_observed':git.stdout.strip() if git.returncode==0 else None,
       'native_jobs_launched':0,'results':{},'errors':{}}
    jobs={'atlas_review':lambda:review(repo,out,'atlas_review',archive_date=archive_date),
          'input_inventory':lambda:inventory(repo,out,'input_inventory',hash_dense=hash_dense),
          'covariance_review':lambda:inspect_covariance(repo,out,'covariance_review'),
          'replication_review':lambda:review_replication(repo,out,'replication_review')}
    if scan_shards:
        from .shards24 import audit_shards
        jobs['shard_audit']=lambda:audit_shards(repo,out,'shard_audit')
    for key,fn in jobs.items():
        try:
            artifact=fn();state['results'][key]=str(artifact)
        except (OSError,ValueError,KeyError,StopIteration) as exc:
            state['errors'][key]=f'{type(exc).__name__}: {exc}'
    rows=[]
    for n,(stage,scope) in STEPS.items():
        status='NOT_EXECUTED_BY_HOST_AUDIT';reason='Bind real production artifacts; this scan is not the analysis itself'
        if n==1:
            status='BLOCKED_SOFTWARE' if env['missing_requirements'] else 'RUNTIME_FOUND_NOT_STATISTICALLY_VALIDATED'
            reason=';'.join(env['missing_requirements']) or 'Run native statistical validation, not just imports'
        elif n==2:
            status='DRAFT_SELECTION_REVIEW_REQUIRED' if 'atlas_review' in state['results'] else 'BLOCKED_ATLAS_INPUT'
            reason='Proposals use measured full-family FDR; no primary pair is silently frozen'
        elif n==3:
            if 'input_inventory' in state['results']:
                r=read_json(out/'input_inventory/summary.json');status='BLOCKED_DENSE_INPUTS' if not r['dense_present'] else 'DENSE_BYTES_REQUIRE_ROW_AND_SCALE_VALIDATION'
                reason=f"{r['dense_present']}/{r['traits']} dense input paths contain actual bytes"
            else:status='BLOCKED_SOURCE_REGISTRY'
        elif n==5:
            status='ARCHIVED_GLOBAL_ARITHMETIC_VERIFIED' if 'atlas_review' in state['results'] else 'BLOCKED_ATLAS_INPUT'
            reason='Existing P/FDR results checked, not a new LDSC calculation'
        elif n==6:
            status='EXISTING_TRACK_B_LEDGER_REVIEWED' if 'replication_review' in state['results'] else 'REPLICATION_LEDGER_UNAVAILABLE'
            reason='Preserved original directional/nonoverlap caveats; no new Brain6 replication run'
        elif n==18:
            status='EXISTING_MODELS_REQUIRE_REVALIDATION' if 'covariance_review' in state['results'] else 'COVARIANCE_EXPORT_UNAVAILABLE'
            reason='Check computed eigenvalues and saved fit flags; no silent smoothing or factor acceptance'
        elif n==24:status='RELEASE_BLOCKED';reason='All required production units must pass the completion24 audit'
        rows.append({'step':n,'stage':stage,'status':status,'reason':reason})
    write_tsv(out/'steps_1_24.tsv',list(rows[0]),rows)
    state['computational_project_complete']=False
    write_json(out/'HOST_AUDIT.json',state)
    return state
