import requests,gzip,struct,re
u='https://zenodo.org/api/records/14908182/files/Ext_eqtl_full_assoc.tsv.gz/content'; size=4595411722
for frac in [0.05,0.2,0.35,0.5,0.65,0.8,0.95]:
 start=int(size*frac);r=requests.get(u,headers={'Range':f'bytes={start}-{start+131071}'},timeout=35);r.raise_for_status();b=r.content;m=re.search(b'\x1f\x8b\x08\x04....\x00\xff\x06\x00BC\x02\x00',b,re.DOTALL)
 if not m: print(frac,'no bgzf');continue
 i=m.start();n=struct.unpack('<H',b[i+16:i+18])[0]+1
 try: data=gzip.decompress(b[i:i+n]).decode();line=data.split('\n')[1]; cols=line.split('\t');print(frac,'start',start+i,'first_complete',cols[0],cols[2],cols[3])
 except Exception as e:print(frac,'ERROR',repr(e))
