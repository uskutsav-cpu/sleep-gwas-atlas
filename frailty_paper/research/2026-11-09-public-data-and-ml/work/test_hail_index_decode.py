import lz4.block,struct
b=open('work/stage6_gnomad_fin_index_part.bin','rb').read();off=0;output=[]
while off<len(b):
 if len(b)-off<4: print('truncated length',len(b)-off);break
 n=struct.unpack_from('<I',b,off)[0];off+=4
 if n==0:print('zero terminator',off);break
 if n>len(b)-off:print('bad length',n,off);break
 raw=struct.unpack_from('<I',b,off)[0];output.append(lz4.block.decompress(b[off+4:off+n],uncompressed_size=raw));off+=n
s=b''.join(output);print('blocks',len(output),'bytes',len(s),'trailing input',len(b)-off,'tail',s[-30:].hex())
i=0;rows=[]
def byte():
 global i
 x=s[i];i+=1;return x
def leb():
 v=0;shift=0
 while 1:
  x=byte();v|=(x&127)<<shift
  if x<128:return v
  shift+=7
  if shift>63:raise ValueError('oversized leb')
def string():
 global i
 n=leb();out=s[i:i+n].decode();i+=n;return out
try:
 while i<len(s):
  marker=byte()
  if marker==0:print('end marker',i);break
  assert marker==1,('marker',marker,i)
  mask=byte();assert mask==0,('rowmask',mask,i)
  chrom=string();pos=leb(); n=leb();mask=byte();assert mask==0 and n<=4,('alleles',n,mask,i)
  al=[string() for _ in range(n)];mask=byte();assert mask==0,('popmask',mask,i)
  ac=leb();af=struct.unpack_from('<d',s,i)[0];i+=8;an=leb();hom=leb();idx=leb(); rows.append((chrom,pos,al,ac,af,an,hom,idx))
 print('parsed',len(rows),'last offset',i,'first',rows[:2],'last',rows[-2:])
except Exception as e:print('ERROR parsed',len(rows),'offset',i,'exception',repr(e),'upcoming',s[max(0,i-15):i+20].hex());raise
