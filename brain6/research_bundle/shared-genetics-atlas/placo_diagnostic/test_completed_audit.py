import tempfile,pathlib,gzip,csv,subprocess,json,math
script=pathlib.Path('outputs/shared-genetics-atlas/placo_diagnostic/audit_completed.py').resolve()
with tempfile.TemporaryDirectory() as tmp:
 root=pathlib.Path(tmp);inp=root/'in';res=root/'results';out=root/'out';inp.mkdir();res.mkdir()
 def write(p,cols,rows):
  with gzip.open(p,'wt') as f:
   w=csv.DictWriter(f,fieldnames=cols,delimiter='\t');w.writeheader();w.writerows(rows)
 cols=['SNP','A1','A2','Z'];s=[dict(SNP='rs1',A1='A',A2='C',Z=1),dict(SNP='rs2',A1='G',A2='T',Z=2),dict(SNP='rs3',A1='A',A2='G',Z=3)];write(inp/'S.sumstats.gz',cols,s)
 rcols=['SNP','A1','A2','Z_insomnia','Z_frailty','P_insomnia','P_frailty','P_PLACO_plus','status']
 for pair in ('FI','FFS'):
  t=[dict(r) for r in s];write(inp/(pair+'.sumstats.gz'),cols,t)
  p=lambda z:math.erfc(abs(z)/math.sqrt(2))
  rows=[dict(SNP=r['SNP'],A1=r['A1'],A2=r['A2'],Z_insomnia=r['Z'],Z_frailty=r['Z'],P_insomnia=p(r['Z']),P_frailty=p(r['Z']),P_PLACO_plus=pp,status=status) for r,pp,status in zip(s,[.001,.1,'NA'],['TESTED','TESTED','NUMERICAL_FAILURE'])]
  write(res/(pair+'_screen.tsv.gz'),rcols,rows)
 control=root/'controls.tsv'
 with control.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=['SNP','A1','A2','Z_insomnia','Z_FI','Z_FFS','P_PLACO_FI','P_PLACO_FFS'],delimiter='\t');w.writeheader();w.writerow(dict(SNP='rs1',A1='A',A2='C',Z_insomnia=1,Z_FI=1,Z_FFS=1,P_PLACO_FI=.001,P_PLACO_FFS=.001))
 cmd=['python3',str(script),'--input',str(inp),'--results',str(res),'--out',str(out),'--expected-n','3','--controls',str(control)];subprocess.run(cmd,check=True,stdout=subprocess.DEVNULL);d=json.loads((out/'completed_audit.json').read_text());assert d['pair_counts']['FI']['threshold_crossing']==1 and d['pair_counts']['FFS']['NUMERICAL_FAILURE']==1 and d['pair_counts']['FFS']['control_matches']==1;assert len(list(csv.DictReader((out/'audited_threshold_crossings.tsv').open(),delimiter='\t')))==2
 # Control values are checked, not merely counted.
 original_control=control.read_text();control.write_text(original_control.replace('0.001','0.002'));p_bad=subprocess.run(cmd,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL);assert p_bad.returncode!=0;control.write_text(original_control)
 # A valid-looking result in a different order must be rejected.
 rows[0]['SNP']='rsWrong';write(res/'FFS_screen.tsv.gz',rcols,rows);p=subprocess.run(cmd,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL);assert p.returncode!=0
 # Missing completed files must not cause partial data to be used.
 (res/'FFS_screen.tsv.gz').unlink();p=subprocess.run(cmd,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL);assert p.returncode!=0
 print('PASS: complete matching, failure retention, family count, control mismatch rejection, SNP mismatch rejection, and missing completion rejection')
