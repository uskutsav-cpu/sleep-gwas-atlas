import sys,zipfile,pathlib,hashlib
from stage5_remotezip_lib import RangeFile
u='https://zenodo.org/api/records/15297965/files/data_MPRA.zip/content';path=sys.argv[1];target=pathlib.Path('/Volumes/Extreme SSD/codex-sleep-frailty-work')/path.split('/')[-1]
with zipfile.ZipFile(RangeFile(u)) as z:data=z.read(path)
target.write_bytes(data);print(path,'bytes',len(data),'sha256',hashlib.sha256(data).hexdigest())
