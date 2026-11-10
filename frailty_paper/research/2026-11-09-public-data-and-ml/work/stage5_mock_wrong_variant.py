import os,runpy,sys,tempfile,pathlib,json
import pandas as pd
from alphagenome.models.v1 import dna_client,variant_scorers
class Meta: pass
class Model:
 def output_metadata(self,**kwargs):return Meta()
 def score_variant(self,**kwargs):return ['mock']
dna_client.create=lambda *args,**kwargs:Model()
variant_scorers.tidy_scores=lambda *args,**kwargs:pd.DataFrame([{'variant_id':'chr1:1:A>C','output_type':'RNA_SEQ','variant_scorer':'synthetic','raw_score':0.1}])
os.environ['STAGE5_TEST_KEY']='SYNTHETIC_NON_SECRET'
with tempfile.TemporaryDirectory() as td:
 store=pathlib.Path(td);sys.argv=['run','--queue','outputs/stage5_reproducible_scoring/score_queue.csv','--limit','1','--key-env','STAGE5_TEST_KEY','--store',td]
 try:runpy.run_path('outputs/stage5_reproducible_scoring/run_alphagenome.py',run_name='__main__');raise AssertionError('incorrectly succeeded')
 except SystemExit as e:assert e.code==1
 receipts=[p for p in store.glob('*.json') if p.name!='run_config.json'];assert len(receipts)==1
 r=json.loads(receipts[0].read_text());assert r['status']=='failed' and 'variant identifiers' in r['error'];assert not list(store.glob('*.csv.gz'))
 print('PASS: mismatched variant rejected before score file is saved')
