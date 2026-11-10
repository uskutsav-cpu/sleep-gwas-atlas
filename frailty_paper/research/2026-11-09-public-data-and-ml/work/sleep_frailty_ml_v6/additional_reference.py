"""Read the official UCSC GRCh38 twoBit file; ignore soft masking only."""
import struct
from common import *
class TwoBit:
    def __init__(self,path):
        self.f=open(path,'rb');sig=self.f.read(4)
        self.order='>' if sig==b'\x1aA\'C' else '<'
        self.f.seek(0);signature,version,n,reserved=self.ints(4)
        assert signature==0x1A412743 and version==0 and reserved==0
        self.index={};self.records={}
        for _ in range(n):
            size=self.f.read(1)[0];name=self.f.read(size).decode();self.index[name]=self.ints(1)[0]
    def ints(self,n):return struct.unpack(self.order+str(n)+'I',self.f.read(n*4)) if n else ()
    def sequence(self,chrom,start,end):
        if chrom not in self.records:
            self.f.seek(self.index[chrom]);size,nn=self.ints(2);ns=self.ints(nn);nl=self.ints(nn);nm=self.ints(1)[0];self.ints(nm);self.ints(nm);assert self.ints(1)[0]==0
            self.records[chrom]=(size,ns,nl,self.f.tell())
        size,ns,nl,offset=self.records[chrom];assert 0<=start<end<=size
        a,b=start//4,(end+3)//4;self.f.seek(offset+a);raw=self.f.read(b-a)
        seq=''.join('TCAG'[(byte>>shift)&3] for byte in raw for shift in [6,4,2,0])[start%4:start%4+end-start]
        seq=list(seq)
        for s,n in zip(ns,nl):
            for i in range(max(start,s),min(end,s+n)):seq[i-start]='N'
        return ''.join(seq)
def get_reference():
    p=SOURCES/'official_hg38.2bit';m=SOURCES/'official_hg38_md5sum.txt'
    if not m.exists():fetch('https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/md5sum.txt',m,100_000)
    if not p.exists():fetch('https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit',p,1_300_000_000)
    expected=next(line.split()[0] for line in m.read_text().splitlines() if line.split()[-1]=='hg38.2bit')
    h=hashlib.md5()
    with p.open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
    assert h.hexdigest()==expected
    genome=TwoBit(p)
    # Independently compare this reader against all existing UCSC JSON queries.
    verified=0
    for cp in (DATA.parent/'sleep_frailty_ml_v5/kosicki_reference_cache').glob('*.json'):
        if cp.name.startswith('._'):continue
        r=json.loads(cp.read_text());assert genome.sequence(r['chrom'],r['start'],r['end'])==r['dna'].upper();verified+=1
    assert verified>=500
    write_json(OUT/'additional_reference_audit.json',{'source_url':'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit','publisher_MD5_matches':True,'sha256':sha(p),'independent_UCSC_JSON_windows_checked':verified,'reader':str(WORK/'additional_reference.py')})
    return genome
if __name__=='__main__':get_reference();print('OFFICIAL REFERENCE CHECKED',flush=True)
