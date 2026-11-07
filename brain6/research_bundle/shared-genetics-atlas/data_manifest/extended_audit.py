import csv,hashlib,json,os
from pathlib import Path
p=Path(__file__).resolve().parent;atlas=p.parent
load=lambda rel:json.loads((atlas/rel).read_text())
records=[]
def add(role,path,expected,source,limitation):records.append({'role':role,'path':str(path),'expected_sha256':expected,'receipt':source,'limitation':limitation})
# Previously audited primary LDSC source and reference inventory.
for r in csv.DictReader((p/'ldsc_inputs.tsv').open(),delimiter='\t'):add(r['role'],r['path'],r['expected_sha256'],'data_manifest/ldsc_inputs.tsv','Processed inputs or code; hash does not authenticate upstream source or complete environment.')
# 32 public NHANES source files.
nh=load('component_exclusion/nhanes_public/phenotype_pilot/input_verification.json')
for x in nh['files']:add('NHANES public XPT',x['path'],x['expected'],'component_exclusion/nhanes_public/phenotype_pilot/input_verification.json','Public source copy; publisher checksum unavailable; provisional phenotype reconstruction.')
raw=load('global_audit/input_n_audit/raw_sample.json')
for t in ('insomnia','adhd'):add('Brain6 raw '+t,raw[t]['raw_path'],raw[t]['raw_sha256'],'global_audit/input_n_audit/raw_sample.json','Saved local hash; upstream provenance and source processing need review.')
base=Path('/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1');screen=load('variant_context/brain6_fi_ffs_screen/receipt.json')
add('Atkins FI raw',base/'data/raw/frailty.txt.gz',screen['hashes']['FI'],'variant_context/brain6_fi_ffs_screen/receipt.json','Current source copy; historical frozen FI bytes differ and unavailable.')
add('Ye FFS raw',base/'frailty_paper/data/gwas/fried_frailty_score/Frailty_GWAS_maf0.01_info0.9_hwe1e-6_FUMA.tsv',screen['hashes']['FFS'],'variant_context/brain6_fi_ffs_screen/receipt.json','Archived source hash; effect-scale and some adapter semantics remain unverified.')
fg=load('replication/finngen_adhd_archived_munged/receipt.json');add('Archived ADHD-derived munged',fg['source'],fg['sha256'],'replication/finngen_adhd_archived_munged/receipt.json','Raw FinnGen origin and ingest receipt unverified; not an independent replication.')
m=load('global_audit/brain6_frailty_covariance/three_frailty_common_support/matching_receipt.json');base2=Path('/Volumes/Brain6Work/brain6-frailty-paired')
for name,expected in m['earlier_input_sha256'].items():add('Three-frailty paired input',name,expected,'global_audit/brain6_frailty_covariance/three_frailty_common_support/matching_receipt.json','Derived common support; ancestry, phenotype and source overlap remain limitations.')
add('Three-frailty HFRS input',base2/'HFRS.sumstats.gz',m['hfrs_input_sha256'],'global_audit/brain6_frailty_covariance/three_frailty_common_support/matching_receipt.json','Derived HFRS support; source N and reference fit remain uncertain.')
x=load('variant_context/rs62264764_brain_xqtl/ldz_readiness/matching_receipt.json')
for key in ('oli','exc'):
 d=x['cells'][key];add('Brain QTL matched '+key,d['file'],d['sha256'],'variant_context/rs62264764_brain_xqtl/ldz_readiness/matching_receipt.json','Derived rows; source data and study-matched LD require independent validation.')
r=load('variant_context/brain6_fi_ffs_screen/three_score_receipt.json');add('UKB HFRS extracted table','/Volumes/Extreme SSD/brain6-hfrs-candidate-results-20261028.tsv',r['hfrs_local_table_sha256'],'variant_context/brain6_fi_ffs_screen/three_score_receipt.json','Prior extraction receipt binds table; raw source not embedded.')
seen=set();summary={}
for x in records:
 assert x['path'] not in seen,x['path'];seen.add(x['path']);q=Path(x['path'])
 if not q.is_file():x.update(bytes='',observed_sha256='',status='missing')
 else:
  try:
   h=hashlib.sha256()
   with q.open('rb') as f:
    for b in iter(lambda:f.read(1048576),b''):h.update(b)
   x.update(bytes=q.stat().st_size,observed_sha256=h.hexdigest(),status='match' if h.hexdigest()==x['expected_sha256'] else 'mismatch')
  except OSError as e:x.update(bytes='',observed_sha256='',status='unreadable:'+type(e).__name__)
 summary[x['status']]=summary.get(x['status'],0)+1
with (p/'external_inputs.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=['role','path','bytes','expected_sha256','observed_sha256','status','receipt','limitation'],delimiter='\t');w.writeheader();w.writerows(records)
(p/'external_check.json').write_text(json.dumps({'counts':summary,'total':len(records),'scope':'Curated key dependencies only; not an exhaustive whole-atlas or scientific provenance audit.','unavailable_or_unverified':['Modified FI excluding sleep and fatigue: no complete public per-SNP file identified','Archived FinnGen ADHD-derived file: raw source and ingest receipt not found','Historical frozen FI bytes differ from current copy and are unavailable','Study-matched LD and complete source metadata for selected molecular and local analyses are not established']},indent=2)+'\n')
print(json.dumps({'counts':summary,'total':len(records)},indent=2))
