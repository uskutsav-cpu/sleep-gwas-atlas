import pathlib,tempfile,csv,gzip,json,hashlib,math,subprocess
root=pathlib.Path('outputs/shared-genetics-atlas/placo_diagnostic').resolve();script=root/'candidate_source_audit.py'
with tempfile.TemporaryDirectory() as tmp:
 d=pathlib.Path(tmp);a=d/'audit';a.mkdir();o=d/'out'
 hits=[{'pair':'FI','SNP':'rs1','A1':'A','A2':'C','Z_insomnia':'2','Z_frailty':'2','p_placo_plus':'1e-9'},{'pair':'FI','SNP':'rs2','A1':'G','A2':'T','Z_insomnia':'2','Z_frailty':'2','p_placo_plus':'2e-9'},{'pair':'FFS','SNP':'rs3','A1':'A','A2':'G','Z_insomnia':'2','Z_frailty':'2','p_placo_plus':'3e-9'}]
 def wr(path,rows,fields=None):
  with (gzip.open(path,'wt') if str(path).endswith('.gz') else open(path,'w')) as f:
   w=csv.DictWriter(f,fieldnames=fields or list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
 p=a/'audited_threshold_crossings.tsv';wr(p,hits);(a/'completed_audit.json').write_text(json.dumps({'hits_sha256':hashlib.sha256(p.read_bytes()).hexdigest()}));pval=math.erfc(2/math.sqrt(2))
 s=[{'SNP':h['SNP'],'A1':h['A1'],'A2':h['A2'],'BETA':.2,'SE':.1,'P':pval} for h in hits];wr(d/'S.tsv',s)
 fi=[{'variant_id':'rs1','effect_allele':'C','other_allele':'A','beta':-.2,'standard_error':.1,'p_value':pval},{'variant_id':'rs2','effect_allele':'G','other_allele':'T','beta':.2,'standard_error':.1,'p_value':pval}];fi.append(fi[-1]);wr(d/'FI.gz',fi)
 wr(d/'FFS.tsv',[{'SNP':'rs3','A1':'A','A2':'T','BETA':.2,'SE':.1,'P':pval}])
 cmd=['python3',str(script),'--audit',str(a),'--s',str(d/'S.tsv'),'--fi',str(d/'FI.gz'),'--ffs',str(d/'FFS.tsv'),'--out',str(o)];subprocess.run(cmd,check=True,stdout=subprocess.DEVNULL)
 r=list(csv.DictReader((o/'candidate_source_matches.tsv').open(),delimiter='\t'));assert [x['sensitivity_eligible'] for x in r]==['YES','NO','NO'];assert r[0]['frailty_orientation']=='-1' and abs(float(r[0]['frailty_z_wald'])-2)<1e-12;assert r[1]['frailty_source_status']=='MULTIPLE_VALID_MATCHES' and r[2]['frailty_source_status']=='NO_MATCH'
 subprocess.run(['Rscript',str(root/'candidate_sensitivity.R'),str(o)],check=True,stdout=subprocess.DEVNULL);s=list(csv.DictReader((o/'candidate_statistic_sensitivity.tsv').open(),delimiter='\t'));assert 0<float(s[0]['p_source_p'])<=1 and s[1]['p_source_p']=='NA'
 (a/'completed_audit.json').write_text(json.dumps({'hits_sha256':'incorrect'}));assert subprocess.run(cmd,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL).returncode!=0
 print('PASS: allele reversal, duplicate and missing exclusion, R calculation, and input hash gate')
