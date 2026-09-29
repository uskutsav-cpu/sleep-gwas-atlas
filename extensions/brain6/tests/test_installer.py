from pathlib import Path
import subprocess
import pytest
from brain6.install import repository_identity,apply_overlay,DEFAULT_BRANCH
from brain6.io import ContractError,write_json,sha256

@pytest.mark.parametrize('remote',[
 'https://github.com/uskutsav-cpu/sleep-gwas-atlas.git',
 'git@github.com:uskutsav-cpu/sleep-gwas-atlas.git',
 'ssh://git@github.com/uskutsav-cpu/sleep-gwas-atlas',
 'https://github.com/uskutsav-cpu/sleep-gwas-atlas'])
def test_remote_identification(remote):
    assert repository_identity(remote)=='uskutsav-cpu/sleep-gwas-atlas'

@pytest.mark.parametrize('remote',[
 'https://github.com.evil.org/uskutsav-cpu/sleep-gwas-atlas.git',
 'http://github.com/uskutsav-cpu/sleep-gwas-atlas',
 'file:///somewhere/repo','https://evil.example/github.com/uskutsav-cpu/sleep-gwas-atlas'])
def test_remote_reject_lookalikes(remote):
    with pytest.raises(ContractError):repository_identity(remote)

@pytest.fixture
def git_fixture(tmp_path):
    repo=tmp_path/'repo';repo.mkdir()
    def run(*args):
        r=subprocess.run(['git','-C',str(repo),*args],capture_output=True,text=True,check=True);return r.stdout.strip()
    run('init','-b','main');run('config','user.name','Synthetic Test Fixture');run('config','user.email','test@example.invalid')
    run('remote','add','origin','https://github.com/uskutsav-cpu/sleep-gwas-atlas.git')
    (repo/'README.md').write_text('original');run('add','README.md');run('commit','-m','test fixture')
    payload=tmp_path/'payload';dest=payload/'extensions/brain6/new.txt';dest.parent.mkdir(parents=True);dest.write_text('addition')
    manifest=tmp_path/'manifest.json';write_json(manifest,{'extensions/brain6/new.txt':sha256(dest)})
    return repo,payload,manifest,run

def test_default_installer_does_not_write(git_fixture):
    r,p,m,git=git_fixture;before=git('rev-parse','HEAD')
    x=apply_overlay(r,p,m)
    assert x['action']=='CHECK_ONLY' and not (r/'extensions').exists()
    assert git('rev-parse','HEAD')==before and git('branch','--show-current')=='main'

def test_installer_commit_preserves_unrelated_changes(git_fixture):
    r,p,m,git=git_fixture;(r/'README.md').write_text('uncommitted user work')
    (r/'unrelated.txt').write_text('untracked user work')
    x=apply_overlay(r,p,m,apply=True,commit=True)
    assert x['commit'] and git('branch','--show-current')==DEFAULT_BRANCH
    assert git('show','--pretty=','--name-only','HEAD')=='extensions/brain6/new.txt'
    assert (r/'README.md').read_text()=='uncommitted user work'
    assert (r/'unrelated.txt').read_text()=='untracked user work'
    assert not x['remote_changed']

def test_installer_refuses_conflict_before_branch(git_fixture):
    r,p,m,git=git_fixture;f=r/'extensions/brain6/new.txt';f.parent.mkdir(parents=True);f.write_text('user file')
    with pytest.raises(ContractError,match='Conflict'):apply_overlay(r,p,m,apply=True)
    assert f.read_text()=='user file' and git('branch','--show-current')=='main'

def test_installer_refuses_staged_work(git_fixture):
    r,p,m,git=git_fixture;(r/'README.md').write_text('staged user work');git('add','README.md')
    with pytest.raises(ContractError,match='staged'):apply_overlay(r,p,m,apply=True,commit=True)
    assert not (r/'extensions').exists()

def test_installer_idempotent_identical_source(git_fixture):
    r,p,m,git=git_fixture
    apply_overlay(r,p,m,apply=True)
    x=apply_overlay(r,p,m,apply=True);assert not x['new_files']

def test_installer_refuses_modified_payload(git_fixture):
    r,p,m,git=git_fixture;(p/'extensions/brain6/new.txt').write_text('corrupted')
    with pytest.raises(ContractError,match='checksum'):apply_overlay(r,p,m)

def test_installer_no_push_without_explicit_commit(git_fixture):
    r,p,m,git=git_fixture
    with pytest.raises(ContractError,match='--push'):apply_overlay(r,p,m,apply=True,push=True)

def test_installer_wrong_repo(git_fixture):
    r,p,m,git=git_fixture;git('remote','set-url','origin','https://github.com/uskutsav-cpu/not-this-project.git')
    with pytest.raises(ContractError,match='Expected'):apply_overlay(r,p,m,apply=True)


def test_exact_prior_payload_upgrade_requires_opt_in(git_fixture,tmp_path):
    r,p,m,git=git_fixture
    old=r/'extensions/brain6/new.txt';old.parent.mkdir(parents=True);old.write_text('known v0.1')
    baseline=tmp_path/'baseline.json';write_json(baseline,{'extensions/brain6/new.txt':sha256(old)})
    with pytest.raises(ContractError,match='Conflict'):apply_overlay(r,p,m,apply=True)
    plan=apply_overlay(r,p,m,baseline_manifest=baseline)
    assert plan['upgrade_files']==['extensions/brain6/new.txt']
    assert old.read_text()=='known v0.1'
    plan=apply_overlay(r,p,m,apply=True,baseline_manifest=baseline)
    assert old.read_text()=='addition' and plan['remote_changed'] is False


def test_upgrade_preserves_local_edits(git_fixture,tmp_path):
    r,p,m,git=git_fixture
    old=r/'extensions/brain6/new.txt';old.parent.mkdir(parents=True);old.write_text('baseline')
    baseline=tmp_path/'baseline.json';write_json(baseline,{'extensions/brain6/new.txt':sha256(old)})
    old.write_text('user changes that must survive')
    with pytest.raises(ContractError,match='Conflict'):apply_overlay(r,p,m,apply=True,baseline_manifest=baseline)
    assert old.read_text()=='user changes that must survive'
    assert git('branch','--show-current')=='main'


def test_upgrade_commit_only_exact_payload(git_fixture,tmp_path):
    r,p,m,git=git_fixture
    old=r/'extensions/brain6/new.txt';old.parent.mkdir(parents=True);old.write_text('old')
    git('add','extensions/brain6/new.txt');git('commit','-m','Synthetic old package')
    baseline=tmp_path/'baseline.json';write_json(baseline,{'extensions/brain6/new.txt':sha256(old)})
    (r/'README.md').write_text('important uncommitted work')
    apply_overlay(r,p,m,apply=True,commit=True,baseline_manifest=baseline)
    assert git('show','--pretty=','--name-only','HEAD')=='extensions/brain6/new.txt'
    assert (r/'README.md').read_text()=='important uncommitted work'


def test_symlink_upgrade_forbidden(git_fixture,tmp_path):
    r,p,m,git=git_fixture
    outside=tmp_path/'elsewhere';outside.write_text('baseline')
    target=r/'extensions/brain6/new.txt';target.parent.mkdir(parents=True);target.symlink_to(outside)
    baseline=tmp_path/'baseline.json';write_json(baseline,{'extensions/brain6/new.txt':sha256(outside)})
    with pytest.raises(ContractError):apply_overlay(r,p,m,apply=True,baseline_manifest=baseline)
    assert outside.read_text()=='baseline'


def test_installer_refuses_unfinished_merge(git_fixture):
    repo,payload,manifest,git=git_fixture
    (repo/'.git/MERGE_HEAD').write_text(git('rev-parse','HEAD')+'\n')
    with pytest.raises(ContractError,match='Finish the existing merge'):
        apply_overlay(repo,payload,manifest,apply=True)
    assert git('branch','--show-current')=='main'
    assert not (repo/'extensions').exists()
