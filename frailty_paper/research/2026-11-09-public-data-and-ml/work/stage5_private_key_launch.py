#!/usr/bin/env python3
"""Receive an API key in a terminal without echo and launch scoring without saving the key."""
import getpass,os,subprocess,sys,pathlib
key=getpass.getpass('AlphaGenome API key (input hidden): ').strip()
if not key:sys.exit('No key entered; no request was made.')
env=os.environ.copy();env['ALPHAGENOME_API_KEY']=key;env['OPENBLAS_NUM_THREADS']='1';env['OMP_NUM_THREADS']='1'
# Start a background process; do not write the key to a file, the command line, or stdout.
log=pathlib.Path('work/stage5_live_scoring.log').open('ab',buffering=0)
child=subprocess.Popen([sys.executable,'work/stage5_run_alphagenome.py'],env=env,stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,start_new_session=True,close_fds=True)
print('AlphaGenome scoring started; process ID:',child.pid,flush=True)
