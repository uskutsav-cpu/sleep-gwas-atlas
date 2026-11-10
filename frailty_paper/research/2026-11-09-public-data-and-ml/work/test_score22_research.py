import csv,subprocess,sys,tempfile,math
from pathlib import Path
base=Path('outputs').resolve();score=base/'score_22_research.py';weights=list(csv.DictReader((base/'unvalidated_22_variant_score_prototype.csv').open())); assert len(weights)==22
fields=['sample_id','variant','dosage_allele','other_allele','dosage']; obs=[]
for i,w in enumerate(weights):
 d=(i%5)*.5
 obs.append(dict(sample_id='direct',variant=w['variant'],dosage_allele=w['counted_allele'],other_allele=w['other_allele'],dosage=d))
 obs.append(dict(sample_id='reversed',variant=w['variant'],dosage_allele=w['other_allele'],other_allele=w['counted_allele'],dosage=2-d))
def write(path,rows):
 with path.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
def run(inp,out,allow_reversal=False):
 args=[sys.executable,str(score),'--dosages',str(inp),'--output',str(out)]
 if allow_reversal: args.append('--allow-verified-reversal')
 return subprocess.run(args,capture_output=True,text=True)
with tempfile.TemporaryDirectory(dir='work') as t:
 p=Path(t);inp=p/'input.csv';out=p/'output.csv';write(inp,obs)
 r=run(inp,out);assert r.returncode!=0 and not out.exists() and 'Reversed dosage' in r.stderr,r.stderr
 r=run(inp,out,allow_reversal=True);assert r.returncode==0,r.stderr
 result={x['sample_id']:x for x in csv.DictReader(out.open())}
 for k,weight in [('FI_raw_sum','FI_beta_weight'),('FFS_raw_sum','FFS_beta_weight_inferred'),('UKB_insomnia_raw_sum','UKB_insomnia_beta_weight')]:
  expected=sum(((i%5)*.5)*float(w[weight]) for i,w in enumerate(weights))
  assert math.isclose(float(result['direct'][k]),expected,abs_tol=1e-14) and math.isclose(float(result['reversed'][k]),expected,abs_tol=1e-14)
 assert run(inp,out,allow_reversal=True).returncode!=0
 mutations={'missing':obs[:-1],'duplicate':obs+[obs[0].copy()],'allele mismatch':[dict(x,other_allele='N') if i==0 else x for i,x in enumerate(obs)],'nonfinite dosage':[dict(x,dosage='nan') if i==0 else x for i,x in enumerate(obs)],'out of range':[dict(x,dosage='2.01') if i==0 else x for i,x in enumerate(obs)]}
 for name,rows in mutations.items():
  write(inp,rows);badout=p/(name.replace(' ','_')+'.csv');r=run(inp,badout,allow_reversal=True);assert r.returncode!=0 and not badout.exists(),(name,r.stdout,r.stderr)
 for name,content in {'duplicate_header':'sample_id,variant,dosage_allele,other_allele,dosage,dosage\na,b,A,G,1,2\n','short_row':'sample_id,variant,dosage_allele,other_allele,dosage\na,b,A,G\n','long_row':'sample_id,variant,dosage_allele,other_allele,dosage\na,b,A,G,1,extra\n'}.items():
  inp.write_text(content);badout=p/(name+'.csv');r=run(inp,badout);assert r.returncode!=0 and not badout.exists() and 'Scoring failed:' in r.stderr,(name,r.stderr)
 print('Passed synthetic direct and reversed allele arithmetic; duplicate, missing, mismatch, nonfinite, range and output-overwrite rejection checks.')
