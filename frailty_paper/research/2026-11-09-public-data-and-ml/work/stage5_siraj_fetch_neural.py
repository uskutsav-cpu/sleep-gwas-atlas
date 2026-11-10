import sys,zipfile,pathlib,hashlib,json,concurrent.futures
from stage5_remotezip_lib import RangeFile
u='https://zenodo.org/api/records/15297965/files/data_MPRA.zip/content';root=pathlib.Path('/Volumes/Extreme SSD/codex-sleep-frailty-work/siraj2026_neural');root.mkdir(exist_ok=True)
with zipfile.ZipFile(RangeFile(u)) as z:members=[x.filename for x in z.infolist() if '/processed_results/DE/' in x.filename and '_SKNSH_' in x.filename and x.filename.endswith('.gz')]
def fetch(name):
 target=root/name.split('/')[-1]
 if not target.exists():
  with zipfile.ZipFile(RangeFile(u)) as z:data=z.read(name)
  target.write_bytes(data)
 data=target.read_bytes();print('downloaded',target.name,len(data),flush=True);return {'member':name,'file':str(target),'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()}
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:r=list(pool.map(fetch,members))
pathlib.Path('outputs/stage5_siraj2026_neural_source_receipt.json').write_text(json.dumps({'url':u,'selection':'all SKNSH processed allele-differential result files','files':r},indent=2)+'\n');print('completed',len(r))
