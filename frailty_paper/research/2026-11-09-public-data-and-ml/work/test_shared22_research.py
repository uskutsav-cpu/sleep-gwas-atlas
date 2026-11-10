import csv,subprocess,sys,tempfile,pathlib,math
p=pathlib.Path('outputs').resolve();model=list(csv.DictReader(open(p/'unvalidated_22_shared_direction_score.csv')));assert len(model)==22
cols=['sample_id','variant','dosage_allele','other_allele','dosage'];obs=[]
for i,r in enumerate(model):
 d=(i%5)*.5
 obs.append(dict(sample_id='direct',variant=r['variant'],dosage_allele=r['counted_allele'],other_allele=r['other_allele'],dosage=d))
 obs.append(dict(sample_id='reversed',variant=r['variant'],dosage_allele=r['other_allele'],other_allele=r['counted_allele'],dosage=2-d))
with tempfile.TemporaryDirectory(dir='work') as d:
 inp=pathlib.Path(d)/'input.csv';out=pathlib.Path(d)/'output.csv'
 with inp.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=cols);w.writeheader();w.writerows(obs)
 cmd=[sys.executable,str(p/'shared_score_22_research.py'),'--dosages',str(inp),'--output',str(out)];r=subprocess.run(cmd,capture_output=True,text=True);assert r.returncode==0,r.stderr
 result=list(csv.DictReader(out.open()));expected=sum((i%5)*.5 for i in range(22));assert len(result)==2 and all(math.isclose(float(x['shared_direction_allele_count']),expected) for x in result)
 assert subprocess.run(cmd,capture_output=True).returncode!=0
 print('Passed direct and reversed allele counts, expected total, and output-overwrite rejection.')
