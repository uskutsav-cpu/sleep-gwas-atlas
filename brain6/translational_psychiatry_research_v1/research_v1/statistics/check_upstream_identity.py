"""Record the public upstream HEAD without cloning or modifying any worktree."""
import datetime, json, subprocess
from calibration import HERE

def main():
    cmd=['git','ls-remote','https://github.com/qlu-lab/SUPERGNOVA.git','HEAD']
    result=subprocess.run(cmd,capture_output=True,text=True,timeout=30)
    record={'checked_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'command':cmd,'exit_code':result.returncode,'stdout':result.stdout,'stderr':result.stderr,'frozen_upstream_commit':'319e84e114a4f954005a4592756c56cfee083667','classification':'PUBLIC_SOURCE_IDENTITY_NOT_EMPIRICAL_VALIDATION'}
    (HERE/'upstream_head_observation.json').write_text(json.dumps(record,indent=2)+'\n')

if __name__=='__main__': main()
