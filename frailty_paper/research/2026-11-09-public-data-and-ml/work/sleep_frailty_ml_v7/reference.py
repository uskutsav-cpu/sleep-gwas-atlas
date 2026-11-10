"""V7 copy of the independently checked twoBit reader; fixed official genome."""
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
    path=V6/"sources/official_hg38.2bit"
    old=json.loads((V6OUT/"additional_reference_audit.json").read_text())
    assert sha(path)==old["sha256"] and old["publisher_MD5_matches"] and old["independent_UCSC_JSON_windows_checked"]>=500
    write_json(OUT/"reference_genome_receipt.json",{"sha256":sha(path),"path":str(path),"previous_reader_source_sha256":sha(ROOT/"work/sleep_frailty_ml_v6/additional_reference.py"),"source_reference_audit_sha256":sha(V6OUT/"additional_reference_audit.json"),"V7_reader_sha256":sha(WORK/"reference.py")})
    return TwoBit(path)
