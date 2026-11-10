import json,os,runpy,sys,tempfile,pathlib,hashlib
import pandas as pd
from alphagenome.models.v1 import dna_client,variant_scorers
runner=pathlib.Path('outputs/stage5_reproducible_scoring/run_alphagenome.py').resolve()
queue=pathlib.Path('outputs/stage5_reproducible_scoring/score_queue.csv').resolve()
state={'calls':0,'empty':False}
class Metadata: pass
class Model:
 def output_metadata(self,**kwargs):return Metadata()
 def score_variant(self,**kwargs):state['calls']+=1;return ['simulated']
dna_client.create=lambda *a,**kw:Model()
def tidy(*a,**kw):
 if state['empty']:return pd.DataFrame()
 return pd.DataFrame([{'variant_id':'chr5:104746478:A>C','output_type':'RNA_SEQ','variant_scorer':'RNA_SEQ','raw_score':0.125}])
variant_scorers.tidy_scores=tidy
os.environ['STAGE5_TEST_KEY']='SYNTHETIC_NON_SECRET'
def run(store):
 sys.argv=[str(runner),'--queue',str(queue),'--limit','1','--key-env','STAGE5_TEST_KEY','--store',str(store)]
 try:runpy.run_path(str(runner),run_name='__main__');return 0
 except SystemExit as e:return e.code or 0
with tempfile.TemporaryDirectory(prefix='stage5-mock-',dir='/Volumes/Extreme SSD/codex-sleep-frailty-work') as td:
 root=pathlib.Path(td);first=root/'success'; assert run(first)==0 and state['calls']==1
 receipts=[p for p in first.glob('*.json') if p.name!='run_config.json' and not p.name.startswith('._')];assert len(receipts)==1
 receipt=json.loads(receipts[0].read_text()); assert receipt['status']=='complete' and receipt['rows']==1 and receipt['returned_output_types']==['RNA_SEQ']
 assert receipt['sha256']==hashlib.sha256((first/receipt['file']).read_bytes()).hexdigest()
 assert run(first)==0 and state['calls']==1
 state['empty']=True;failed=root/'empty';assert run(failed)==1 and state['calls']==2
 receipts=[p for p in failed.glob('*.json') if p.name!='run_config.json' and not p.name.startswith('._')];assert len(receipts)==1
 r=json.loads(receipts[0].read_text());assert r['status']=='failed' and 'no tidy data' in r['error'];assert not [p for p in failed.glob('*.csv.gz') if not p.name.startswith('._')]
 print('MOCK CHECK PASS: success receipt and hash; resume avoids duplicate request; empty response fails without score file')
