import sys,zipfile,pathlib,hashlib
from stage5_remotezip_lib import RangeFile
u='https://zenodo.org/api/records/15297965/files/data_preprocess.zip/content';name=sys.argv[1];target=pathlib.Path('work')/name
with zipfile.ZipFile(RangeFile(u)) as z:data=z.read('data/preprocess/'+name)
target.write_bytes(data);print(name,len(data),hashlib.sha256(data).hexdigest())
