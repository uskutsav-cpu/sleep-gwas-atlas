import pathlib,subprocess,concurrent.futures,hashlib,os
base='https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/'
files={'PGS005228':(25712794,'584c4839ba4d8caa299ff558ed161259'),'PGS005229':(28779650,'168e5c11469c0b27668a4c4d733276b1')}
partdir=pathlib.Path('work/pgs_parts');partdir.mkdir(exist_ok=True);chunk=2097152
jobs=[]
for name,(size,md5) in files.items():
 for start in range(0,size,chunk):jobs.append((name,start,min(start+chunk-1,size-1)))
def fetch(j):
 name,start,end=j;dest=partdir/(name+'_'+str(start)+'.part');url=base+name+'/ScoringFiles/'+name+'.txt.gz?part='+str(start)
 for attempt in range(2):
  result=subprocess.run(['curl','--fail','--silent','--show-error','--max-time','120','--range',str(start)+'-'+str(end),'-o',str(dest),url],capture_output=True)
  if result.returncode==0 and dest.stat().st_size==end-start+1:return name,start,'ok'
 return name,start,'FAILED '+str(result.stderr[:150])
with concurrent.futures.ThreadPoolExecutor(max_workers=12) as e:
 for n,s,state in e.map(fetch,jobs):print(n,s,state,flush=True)
for name,(size,expected) in files.items():
 out=pathlib.Path('work')/(name+'_complete.txt.gz');h=hashlib.md5()
 with out.open('wb') as w:
  for start in range(0,size,chunk):
   data=(partdir/(name+'_'+str(start)+'.part')).read_bytes();w.write(data);h.update(data)
 assert out.stat().st_size==size and h.hexdigest()==expected,(name,out.stat().st_size,h.hexdigest())
 print('VERIFIED',name,'bytes',size,'md5',h.hexdigest(),flush=True)
