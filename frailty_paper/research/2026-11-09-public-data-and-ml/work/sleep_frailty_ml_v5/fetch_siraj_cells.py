"""Retrieve only additional fixed-library cell readouts from the public ZIP."""
import io,zipfile,concurrent.futures
from common import *
URL='https://zenodo.org/api/records/15297965/files/data_MPRA.zip/content'
class RangeFile(io.RawIOBase):
    def __init__(self):
        self.pos=0
        with self.request(0,0) as r:
            assert r.status==206;self.size=int(r.headers['Content-Range'].split('/')[-1]);r.read()
    def request(self,start,end):
        return urllib.request.urlopen(urllib.request.Request(URL,headers={'Range':f'bytes={start}-{end}','User-Agent':'SleepFrailtyResearch/1.0'}),context=SSL,timeout=90)
    def readable(self):return True
    def seekable(self):return True
    def tell(self):return self.pos
    def seek(self,n,whence=0):
        self.pos=({0:0,1:self.pos,2:self.size}[whence]+n);return self.pos
    def read(self,n=-1):
        if n<0:n=self.size-self.pos
        if n==0 or self.pos>=self.size:return b''
        end=min(self.pos+n,self.size)-1
        with self.request(self.pos,end) as r:
            assert r.status==206;data=r.read();assert len(data)==end-self.pos+1
        self.pos=end+1;return data
def main():
    with zipfile.ZipFile(RangeFile()) as z:
        names=[r.filename for r in z.infolist() if '/processed_results/DE/' in r.filename and '/OL42_' in r.filename and '_SKNSH_' not in r.filename and r.filename.endswith('.gz')]
    assert names;print('Additional OL42 members',names,flush=True)
    def one(name):
        target=SOURCES/name.split('/')[-1]
        if not target.exists():
            with zipfile.ZipFile(RangeFile()) as z:b=z.read(name)
            tmp=target.with_suffix('.tmp');tmp.write_bytes(b);tmp.replace(target)
        print('Retrieved',target.name,target.stat().st_size,flush=True)
        return {'member':name,'file':str(target),'bytes':target.stat().st_size,'sha256':sha(target),'ZIP_CRC_checked':True}
    with concurrent.futures.ThreadPoolExecutor(3) as pool:rows=list(pool.map(one,names))
    write_json(SOURCES/'siraj_additional_cell_receipt.json',{'url':URL,'selection':'All OL42 non-SKNSH processed DE files, chosen before viewing their performance','files':rows})
if __name__=='__main__':main()
