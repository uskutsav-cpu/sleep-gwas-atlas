#!/usr/bin/env python3
"""Stream a stored gzip ZIP member in checked HTTP ranges, retrying each range."""
import subprocess,sys
url='https://humandbs.dbcls.jp/files/hum0197/hum0197.v3.BBJ.Ins.v1.zip'
start=109
end=start+1522589138
chunk=16*1024*1024
for lo in range(start,end,chunk):
 hi=min(lo+chunk,end)-1
 for attempt in range(4):
  r=subprocess.run(['curl','-fsSL','--max-time','90','--range',f'{lo}-{hi}',url],stdout=subprocess.PIPE)
  if r.returncode==0 and len(r.stdout)==hi-lo+1:
   sys.stdout.buffer.write(r.stdout);sys.stdout.buffer.flush();break
  print(f'range {lo}-{hi}: retry {attempt+1}',file=sys.stderr,flush=True)
 else:
  print('Range could not be retrieved completely',file=sys.stderr);sys.exit(1)
 if (lo-start)//chunk%20==0: print(f'Completed {(hi-start+1)/1522589138:.0%} of compressed member',file=sys.stderr,flush=True)
