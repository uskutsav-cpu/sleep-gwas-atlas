import csv,gzip,io,json,subprocess,hashlib
from pathlib import Path
root=Path('.')
def read(p):
 with open(p,newline='') as f:return list(csv.DictReader(f))
a={r['marker']:(int(r['chr38']),int(r['pos38'])) for r in read('outputs/stage5_marker_grch37_grch38_liftover_audit.csv') if r['status']=='unique_expected_positive_strand'}
for r in read('outputs/stage6_tpmi_sleep_candidate_inventory.csv'):
 if r['ensembl_liftover_check']=='unique_positive_match':a[r['variant']]=(int(float(r['chr38'])),int(float(r['pos38'])))
for v in ('rs113851554','rs78206187'):
 m=json.loads(Path(f'work/stage6_ensembl_{v}_mapping.json').read_text())['mappings'];assert len(m)==1 and m[0]['mapped']['strand']==1
 a[v]=(int(m[0]['mapped']['seq_region_name']),int(m[0]['mapped']['start']))
files={'score15':'outputs/unvalidated_15_variant_score_prototype.csv','count22':'outputs/unvalidated_22_shared_direction_score.csv'}; candidates={}
for name,f in files.items():
 for r in read(f):
  c=candidates.setdefault(r['variant'],{'variant':r['variant'],'chromosome':a[r['variant']][0],'position38_inferred':a[r['variant']][1],'sets':[],'alleles':set()})
  c['sets'].append(name);c['alleles'].add(frozenset([r['counted_allele'],r['other_allele']]))
assert len(candidates)==33
bypos={}
for c in candidates.values():bypos.setdefault((c['chromosome'],c['position38_inferred']),[]).append(c)
url='https://my.locuszoom.org/gwas/922558/data/?token=fbe996cb953248e49fa54acd9f400da1'
proc=subprocess.Popen(['curl','-fsSL','--max-time','240',url],stdout=subprocess.PIPE)
out=[];n=0
with gzip.GzipFile(fileobj=proc.stdout) as z:
 reader=csv.DictReader(io.TextIOWrapper(z),delimiter='\t');assert reader.fieldnames==['#chrom','pos','rsid','ref','alt','neg_log_pvalue','beta','stderr_beta','alt_allele_freq']
 for r in reader:
  n+=1
  for c in bypos.get((int(r['#chrom']),int(r['pos'])),[]):
   pair=frozenset([r['ref'],r['alt']]);out.append({'candidate':c['variant'],'score_sets':'|'.join(c['sets']),'candidate_pair_match':pair in c['alleles'],**r})
code=proc.wait();assert code==0
assert any(r['candidate']=='rs113851554' and r['ref']=='G' and r['alt']=='T' and abs(float(r['beta'])-.101)<0.00011 for r in out)
p=Path('outputs/stage6_kunorozva_meta_candidate_extract.csv')
with p.open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=['candidate','score_sets','candidate_pair_match','#chrom','pos','rsid','ref','alt','neg_log_pvalue','beta','stderr_beta','alt_allele_freq']);w.writeheader();w.writerows(out)
matched={r['candidate'] for r in out if r['candidate_pair_match']};receipt={'url':url,'source_expected_content_bytes':301671548,'rows_scanned':n,'position_rows_found':len(out),'unique_pair_matches':len(matched),'unmatched_candidates':sorted(set(candidates)-matched),'output_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'limitations':'Discovery sample overlap; study-specific beta allele convention not independently verified; exact allele pair and coordinate matching is not individual score validation.'}
Path('outputs/stage6_kunorozva_meta_candidate_extract_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
