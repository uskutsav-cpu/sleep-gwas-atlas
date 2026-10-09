#!/usr/bin/env python3
"""Complete only the explicitly approved seven-file relocation, retaining path access."""
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess

ROOT=Path(__file__).resolve().parents[2]
P=ROOT/'sleep_unified_research_v4'
ORIGINAL_REPO=Path('/Users/swethasunilkumar/Documents/Codex/2026-09-19/can/work/sleep-gwas-atlas-v03')

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(4*1024**2),b''):h.update(block)
    return h.hexdigest()

def status():
    return subprocess.run(['git','status','--porcelain=v1','-z'],cwd=ORIGINAL_REPO,capture_output=True,check=True).stdout

def main():
    proof_path=P/'logs/internal_archive_relocation_candidate_receipt_v4.json'
    proof=json.loads(proof_path.read_text());plan=json.loads((P/'manifests/internal_archive_relocation_candidate_plan_v4.json').read_text())
    source=Path(plan['source_directory']);target=Path(plan['target_directory'])
    assert source==ORIGINAL_REPO/'data/raw.local-preserved/.archives'
    assert not source.is_symlink() and source.is_dir() and target.is_dir()
    names={r['name'] for r in proof['files']}
    assert {f.name for f in source.iterdir()}==names and len(names)==7
    tracked=subprocess.run(['git','ls-files','--','data/raw.local-preserved/.archives'],cwd=ORIGINAL_REPO,capture_output=True,check=True).stdout
    assert not tracked, 'TRACKED_SOURCE_CANNOT_RELOCATE'
    handles=subprocess.run(['lsof','+D',str(source)],capture_output=True,text=True)
    assert handles.returncode==1 and not handles.stdout.strip(), 'OPEN_SOURCE_HANDLES_OR_LSOF_FAILURE'
    assert not handles.stderr.strip(), 'HANDLE_INSPECTION_ERROR'
    checks=[]
    for row in proof['files']:
        old,new=source/row['name'],target/row['name']
        assert old.is_file() and not old.is_symlink() and new.is_file()
        a,b=sha(old),sha(new)
        assert a==b==row['expected_sha256']
        st=old.stat()
        checks.append({**row,'current_source_sha256':a,'current_target_sha256':b,
                       'original_mtime_ns':st.st_mtime_ns,'original_mode':st.st_mode})
    before=status();free_before=shutil.disk_usage('/System/Volumes/Data').free
    backup=source.with_name('.archives.relocation-pending-v4')
    assert not backup.exists()
    source.rename(backup)
    removed=[]
    try:
        source.symlink_to(target,target_is_directory=True)
        assert source.resolve()==target.resolve()
        for row in checks:
            assert sha(source/row['name'])==row['expected_sha256']
        assert status()==before, 'GIT_WORKTREE_STATUS_CHANGED'
        assert {f.name for f in backup.iterdir()}==names
        for row in checks:
            (backup/row['name']).unlink();removed.append(row['name'])
        backup.rmdir()
    except BaseException:
        if source.is_symlink():source.unlink()
        for name in removed:shutil.copy2(target/name,backup/name)
        if backup.exists():backup.rename(source)
        raise
    after=status();free_after=shutil.disk_usage('/System/Volumes/Data').free
    result={'completed_utc':datetime.now(timezone.utc).isoformat(),
            'user_authorization_exact':'Approve this exact relocation',
            'scope':'Only the seven verified files in data/raw.local-preserved/.archives; no other research project changes authorized.',
            'candidate_receipt_sha256':sha(proof_path),'script_sha256':sha(Path(__file__)),
            'source_alias':str(source),'ssd_target':str(target),'source_alias_resolves_to_target':source.resolve()==target.resolve(),
            'files':checks,'files_relocated':len(checks),'logical_bytes_relocated':sum(r['bytes'] for r in checks),
            'internal_free_bytes_before':free_before,'internal_free_bytes_after':free_after,
            'observed_free_change_bytes':free_after-free_before,'all_original_bytes_preserved':True,
            'git_status_sha256_before':hashlib.sha256(before).hexdigest(),
            'git_status_sha256_after':hashlib.sha256(after).hexdigest(),'git_worktree_status_unchanged':before==after,
            'open_source_handles_before':False,'full_native_internal_floor_bytes':3*1024**3,
            'full_native_internal_gate_pass':free_after>=3*1024**3,
            'metadata_qualification':'Original byte hashes and original metadata are recorded; exFAT does not preserve all APFS inode/permission attributes.',
            'rollback':'After confirming capacity, copy these seven SSD files into a new internal sibling directory, verify all hashes, and atomically replace only the symlink with that directory. Keep SSD originals.'}
    out=P/'logs/authorized_archive_relocation_receipt_v4.json';out.write_text(json.dumps(result,indent=2)+'\n')
    (target.parent/'authorized_archive_relocation_receipt_v4.json').write_bytes(out.read_bytes())
    print(json.dumps({k:result[k] for k in ['files_relocated','internal_free_bytes_before','internal_free_bytes_after','git_worktree_status_unchanged','full_native_internal_gate_pass']}),flush=True)

if __name__=='__main__':main()
