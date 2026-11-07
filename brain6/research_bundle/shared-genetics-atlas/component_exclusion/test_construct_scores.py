import tempfile,pathlib,json,csv,subprocess,math
script=pathlib.Path(__file__).with_name('construct_scores.py')
with tempfile.TemporaryDirectory() as tmp:
 p=pathlib.Path(tmp);fi=['sleep','fatigue']+[f'd{i:02}' for i in range(47)];ffs=['exhaustion','weakness','walking','activity','weight'];m={'fi_deficits':fi,'fi_sleep':'sleep','fi_fatigue':'fatigue','ffs_components':ffs,'ffs_exhaustion':'exhaustion','id_columns':['id']};(p/'mapping.json').write_text(json.dumps(m))
 fields=['id']+fi+ffs;row={k:'0' for k in fields};row.update(id='a',sleep='1',fatigue='0.5',exhaustion='1',walking='1');missing=dict(row,id='b',d00='NA')
 def write(rows):
  with (p/'input.tsv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)
 write([row,missing]);cmd=['python3',str(script),'--input',str(p/'input.tsv'),'--mapping',str(p/'mapping.json'),'--output',str(p/'out.tsv')];q=subprocess.run(cmd,capture_output=True,text=True,check=True);counts=json.loads(q.stdout);assert counts=={'complete':1,'input':2,'missing':1};result=list(csv.DictReader((p/'out.tsv').open(),delimiter='\t'))[0]
 assert math.isclose(float(result['FI49_raw']),1.5/49) and math.isclose(float(result['FI48_no_sleep_raw']),.5/48) and float(result['FI47_no_sleep_fatigue_raw'])==0 and float(result['FFS5_raw'])==2 and float(result['FFS4_no_exhaustion_raw'])==1
 assert subprocess.run(cmd,capture_output=True).returncode!=0,'Existing output must not be overwritten';(p/'out.tsv').unlink();write([row,dict(row,id='c',weakness='0.5')]);assert subprocess.run(cmd,capture_output=True).returncode!=0 and not (p/'out.tsv').exists();write([row,row]);assert subprocess.run(cmd,capture_output=True).returncode!=0 and not (p/'out.tsv').exists()
 print('PASS: fixed complete cases, denominators, exhaustion removal, no overwrite, invalid and duplicate rejection')
