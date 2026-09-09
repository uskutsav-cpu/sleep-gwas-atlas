"""Non-destructive, explicit Git installer for an additive source overlay."""
from __future__ import annotations
import argparse
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
from .io import ContractError,require,sha256,safe_write_path,read_json

EXPECTED_REPOSITORY='uskutsav-cpu/sleep-gwas-atlas'
DEFAULT_BRANCH='feature/brain6-deep-analysis-v02'


def git(repo,*args):
    result=subprocess.run(['git','-C',str(repo),*args],capture_output=True,text=True,check=False)
    require(result.returncode==0,f'Git command failed ({args[0]}): {result.stderr.strip()}')
    return result.stdout.strip()


def repository_identity(remote):
    # Accept HTTPS and SSH GitHub remotes, but not substring matches or lookalike hosts.
    match=re.fullmatch(r'(?:https://(?:[^/@]+@)?github\.com/|git@github\.com:|ssh://git@github\.com/)([^/]+/[^/]+?)(?:\.git)?/?',remote)
    require(match is not None,'Origin must be a GitHub HTTPS/SSH repository URL')
    return match[1]


def inspect_overlay(repo,payload,manifest,*,baseline_manifest=None):
    repo=Path(repo).resolve();payload=Path(payload).resolve()
    require(shutil.which('git') is not None,'Git is required')
    root=Path(git(repo,'rev-parse','--show-toplevel')).resolve()
    require(root==repo,'Pass the repository root, not a subdirectory')
    identity=repository_identity(git(repo,'remote','get-url','origin'))
    require(identity==EXPECTED_REPOSITORY,f'Expected {EXPECTED_REPOSITORY}, not {identity}')
    inventory=read_json(manifest)
    baseline=read_json(baseline_manifest) if baseline_manifest else {}
    actual={str(p.relative_to(payload)) for p in payload.rglob('*') if p.is_file()}
    require(actual==set(inventory),'Overlay file list differs from the delivery checksum inventory')
    changes=[];identical=[];upgrades=[]
    for relative,digest in sorted(inventory.items()):
        source=payload/relative
        require(not source.is_symlink() and all(not p.is_symlink() for p in source.parents if p!=payload.parent),'Payload symlinks are forbidden')
        require(sha256(source)==digest,f'Payload checksum mismatch: {relative}')
        require(relative.startswith('extensions/brain6/') or relative=='.github/workflows/brain6-ci.yml',
                f'Path outside allowed additive extension: {relative}')
        logical=repo/relative
        require(all(not parent.is_symlink() for parent in logical.parents if parent!=repo.parent),
                'Target directory symlinks are not writable by installer')
        target=safe_write_path(repo,relative)
        require(not target.is_symlink(),'Target symlink is not writable by installer')
        if target.exists():
            require(target.is_file(),f'Conflict: {relative} is not a file')
            existing=sha256(target)
            if existing==digest:identical.append(relative)
            else:
                require(relative in baseline and existing==baseline[relative],
                        f'Conflict: {relative}; not an exact recognized baseline, no overwrite')
                upgrades.append(relative)
        else:changes.append(relative)
    return {'repo':str(repo),'repository':identity,'new_files':changes,'upgrade_files':upgrades,'identical_files':identical,'payload_files':sorted(inventory)}


def apply_overlay(repo,payload,manifest,*,apply=False,branch=DEFAULT_BRANCH,commit=False,push=False,baseline_manifest=None):
    require(not commit or apply,'--commit requires --apply')
    require(not push or commit,'--push requires --commit and --apply')
    require(re.fullmatch(r'feature/brain6[-a-zA-Z0-9_/]*',branch) is not None,'Use a dedicated feature/brain6 branch')
    plan=inspect_overlay(repo,payload,manifest,baseline_manifest=baseline_manifest)
    if not apply:return {**plan,'action':'CHECK_ONLY','remote_changed':False}
    repo=Path(repo).resolve();payload=Path(payload).resolve()
    # Existing staged changes must not be committed together with our additions.
    if commit:require(not git(repo,'diff','--cached','--name-only'),'Existing staged changes: commit/unstage those separately; nothing changed')
    current=git(repo,'branch','--show-current')
    require(current,'Detached HEAD: explicitly choose a working branch first')
    if current!=branch:
        branches=git(repo,'for-each-ref','--format=%(refname:short)','refs/heads').splitlines()
        require(branch not in branches,'Target branch already exists; switch to it explicitly and rerun')
        git(repo,'switch','-c',branch)
    baseline=read_json(baseline_manifest) if baseline_manifest else {}
    for relative in plan['new_files']+plan['upgrade_files']:
        source=payload/relative;target=safe_write_path(repo,relative)
        target.parent.mkdir(parents=True,exist_ok=True)
        # Atomic create without replacement; concurrent file creation is a conflict.
        fd,temp=tempfile.mkstemp(prefix='.brain6-install-',dir=target.parent)
        try:
            with os.fdopen(fd,'wb') as f,source.open('rb') as g:
                shutil.copyfileobj(g,f,1024*1024);f.flush();os.fsync(f.fileno())
            if relative in plan['upgrade_files']:
                # Recheck immediately before atomic replacement. Only our known
                # prior payload is replaceable; never a user's edited file.
                require(target.is_file() and not target.is_symlink() and sha256(target)==baseline[relative],
                        f'Concurrent target modification: {relative}; upgrade stopped')
                os.replace(temp,target)
            else:os.link(temp,target)
        finally:
            Path(temp).unlink(missing_ok=True)
    commit_sha=None
    if commit:
        # Use exact paths, never `git add .`, reset, stash, or rewrite someone else's history.
        for i in range(0,len(plan['payload_files']),50):
            git(repo,'add','--',*plan['payload_files'][i:i+50])
        staged=set(git(repo,'diff','--cached','--name-only').splitlines())
        require(staged<=set(plan['payload_files']),'Unrelated concurrent staging detected; refusing commit')
        if staged:
            git(repo,'diff','--cached','--check')
            git(repo,'commit','-m','Add audited Brain6 v0.2 analysis and integrity checks')
        commit_sha=git(repo,'rev-parse','HEAD')
    if push:git(repo,'push','-u','origin',f'HEAD:refs/heads/{branch}')
    return {**plan,'action':'APPLIED','branch':branch,'commit':commit_sha,'remote_changed':push}


def main(argv=None,*,payload_root=None,manifest_path=None,baseline_path=None):
    p=argparse.ArgumentParser(description='Safely add Brain6; defaults to read-only conflict checks.')
    p.add_argument('repo',help='Existing sleep-gwas-atlas checkout root')
    p.add_argument('--apply',action='store_true');p.add_argument('--commit',action='store_true')
    p.add_argument('--upgrade',action='store_true',help='Allow replacement only of exact recognized v0.1 payload files')
    p.add_argument('--push',action='store_true');p.add_argument('--branch',default=DEFAULT_BRANCH)
    a=p.parse_args(argv)
    try:
        require(payload_root is not None and manifest_path is not None,'Run apply_brain6.py from the delivery folder')
        result=apply_overlay(a.repo,payload_root,manifest_path,apply=a.apply,branch=a.branch,commit=a.commit,push=a.push,baseline_manifest=baseline_path if a.upgrade else None)
        print(json.dumps(result,indent=2));return 0
    except (ContractError,OSError,ValueError) as e:
        print(f'brain6 installer: {e}',file=sys.stderr);return 2
