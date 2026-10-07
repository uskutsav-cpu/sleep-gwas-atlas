#!/usr/bin/env python3
"""Frozen clinical-phenotype sensitivity; no variant selection or dense-region claim.

Freeze before fetch. Public browser API only, no bulk cloud files/access forms.
Missing results retain planned slots and are never converted to P=1 or zero effect.
"""
import argparse,csv,hashlib,json,math,subprocess,time
from pathlib import Path

ROOT=Path(__file__).resolve().parents[4]
AREA=ROOT/'brain6/translational_psychiatry_research_v1/research_v1'
DEST=AREA/'sources'
RAW=ROOT/'work/tp-finngen-fixed-20261007'
META=ROOT/'work/tp-literature-20261007'
VARIANTS={
 'A':dict(rsid='rs7105462',chrom='11',pos37=112912048,pos38=113041326,ref='G',alt='A',source_a1='G',source_a2='A'),
 'B':dict(rsid='rs9485410',chrom='6',pos37=101261007,pos38=100813131,ref='T',alt='C',source_a1='T',source_a2='C'),
 'C':dict(rsid='rs77960',chrom='5',pos37=103964585,pos38=104628884,ref='G',alt='A',source_a1='A',source_a2='G'),
}
SLOTS=[('A_INS_ADHD','A','F5_INSOMNIA','F5_ADHD'),('A_INS_MDD','A','F5_INSOMNIA','F5_DEPRESSIO'),('B_INS_MDD','B','F5_INSOMNIA','F5_DEPRESSIO'),('C_INS_ADHD','C','F5_INSOMNIA','F5_ADHD')]
ALPHA=.05/4
MAX_BYTES=5_000_000
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dump(p,d):p.write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')
def tsv(p,rows,fields=None):
 if fields is None:fields=list(rows[0])
 with p.open('w',newline='') as f:
  w=csv.DictWriter(f,fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(rows)
def timestamp():return time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())
def freeze():
 target=DEST/'FIXED_VARIANT_FINGEN_MANIFEST.json'
 if target.exists():raise SystemExit('Manifest already frozen; refusing overwrite')
 rows=[];hashes={}
 for region,v in VARIANTS.items():
  pair={v['ref'],v['alt']}
  for build,pos in [('GRCh37',v['pos37']),('GRCh38',v['pos38'])]:
   p=META/(v['rsid']+'_'+build.lower()+'.json')
   if not p.exists():p=ROOT/'work/literature-20261007'/p.name
   d=json.loads(p.read_text())
   m=[m for m in d['mappings'] if m['assembly_name']==build and m['seq_region_name']==v['chrom'] and m['start']==pos]
   assert len(m)==1 and m[0]['strand']==1 and pair<=set(m[0]['allele_string'].split('/'))
   if build=='GRCh38':assert m[0]['allele_string'].split('/')[0]==v['ref']
   hashes[p.name]=digest(p)
  key=f"{v['chrom']}-{v['pos38']}-{v['ref']}-{v['alt']}"
  rows.append(dict(region=region,**v,key=key,url='https://r13.finngen.fi/api/variant/'+key,max_bytes=MAX_BYTES))
 for n in ['finngen_current_access.html','finngen_r13_config.js','finngen_r13_phenotypes.json','finngen_public_server.py','finngen_public_db.py','finngen_gwas.html','finngen_data_download.html']:
  hashes[n]=digest(META/n)
 for n in ['FIXED_VARIANT_FINGEN_PROTOCOL.md','fixed_variant_finngen.py']:
  hashes[n]=digest(DEST/n)
 phenos={p['phenocode']:p for p in json.loads((META/'finngen_r13_phenotypes.json').read_text())}
 chosen={code:{k:phenos[code].get(k) for k in ['phenocode','phenostring','num_cases','num_controls']} for code in ['F5_INSOMNIA','F5_ADHD','F5_DEPRESSIO']}
 manifest=dict(frozen_utc=timestamp(),run_id='clinical_fixed_variant_finngen_r13_20261007_v1',access='ORDINARY_ADVERTISED_PUBLIC_BROWSER_API',outcomes_seen_before_freeze=False,source_release='FinnGen R13',genome_build='GRCh38',effect_allele='ALT',variants=rows,slots=[dict(slot_id=s,region=r,sleep_code=a,outcome_code=b) for s,r,a,b in SLOTS],family_size=4,alpha_family=.05,alpha_per_slot=ALPHA,test='intersection-union: max(native_P_sleep,native_P_outcome)',missing_rule='NOT_ESTIMATED; retain all four planned slots',classification='PHENOTYPE_SENSITIVITY; discovery_cohort_independence_UNVERIFIED',source_metadata=chosen,source_and_code_sha256=hashes,discovery_allele_evidence=dict(path='work/research-discovery-raw/insomnia.txt.gz',sha256='32848cd92a6324c9048cad4a804cf1546299af1d0281981c6ca43c2046a9065b',receipt='discovery/acquisition_receipts.json',verified_by='Parent agent prior to freeze; original source UNIQUE_ID/A1/A2, not an outcome-directed ALT choice'))
 dump(target,manifest)
 print('FROZEN',digest(target),manifest['frozen_utc'])
def acquire():
 p=DEST/'FIXED_VARIANT_FINGEN_MANIFEST.json';m=json.loads(p.read_text())
 assert digest(DEST/'fixed_variant_finngen.py')==m['source_and_code_sha256']['fixed_variant_finngen.py']
 RAW.mkdir(parents=True,exist_ok=True);receipts=[]
 for v in m['variants']:
  dst=RAW/(v['rsid']+'.json');rec=dict(rsid=v['rsid'],key=v['key'],url=v['url'],manifest_sha256=digest(p),retrieved_utc=timestamp(),max_bytes=v['max_bytes'],raw_path=str(dst.relative_to(ROOT)))
  if dst.exists():raise SystemExit('Refusing silently reuse outcomes without receipt review')
  tmp=dst.with_suffix('.download')
  c=subprocess.run(['curl','--location','--silent','--show-error','--max-time','40','--max-filesize',str(v['max_bytes']),'--output',str(tmp),'--write-out','%{json}',v['url']],capture_output=True,text=True)
  try:meta=json.loads(c.stdout)
  except ValueError:meta={}
  rec.update(http_status=meta.get('http_code'),content_type=meta.get('content_type'),final_url=meta.get('url_effective'),curl_exit=c.returncode)
  rec['bytes']=tmp.stat().st_size if tmp.exists() else 0
  if c.returncode or meta.get('http_code')!=200 or rec['bytes']>v['max_bytes']:
   rec.update(status='UNAVAILABLE',error=c.stderr.strip());tmp.unlink(missing_ok=True)
  else:
   try:d=json.loads(tmp.read_text());assert isinstance(d,dict)
   except (ValueError,AssertionError):rec['status']='INVALID_RESPONSE';tmp.unlink()
   else:tmp.rename(dst);rec.update(status='RETRIEVED',sha256=digest(dst))
  receipts.append(rec);dump(DEST/'fixed_variant_acquisition_receipts.json',receipts)
  print(v['rsid'],rec['status'],rec['bytes'],flush=True)
 analyze()
def finite(x):return isinstance(x,(int,float)) and not isinstance(x,bool) and math.isfinite(x)
def analyze():
 m=json.loads((DEST/'FIXED_VARIANT_FINGEN_MANIFEST.json').read_text());native=[];by={};qcs=[]
 needed={r:{a for _,rr,s,o in SLOTS if rr==r for a in [s,o]} for r in VARIANTS}
 for v in m['variants']:
  file=RAW/(v['rsid']+'.json');d=json.loads(file.read_text()) if file.exists() else {}
  vv=d.get('variant',{});identity=all(str(vv.get(k))==str(v[x]) for k,x in [('chrom','chrom'),('pos','pos38'),('ref','ref'),('alt','alt')])
  # rsid list may be comma-delimited; assembly+alleles uniquely define our point.
  qc=dict(region=v['region'],rsid=v['rsid'],key=v['key'],identity_pass=identity,response_rsids=vv.get('rsids'),result_count=len(d.get('results',[])),response_keys=';'.join(sorted(d)),source_sha256=digest(file) if file.exists() else '',coverage='EXACT_POINT_ONLY; dense_region_completeness_UNKNOWN')
  qcs.append(qc)
  for code in sorted(needed[v['region']]):
   rr=[x for x in d.get('results',[]) if x.get('phenocode')==code]
   z=rr[0] if len(rr)==1 and identity else {}
   reason='ESTIMATED' if z else ('IDENTITY_FAILED' if d and not identity else 'MISSING_EXACT_POINT_OR_PHENOTYPE')
   required=all(finite(z.get(k)) for k in ['pval','beta','sebeta','maf'])
   required=required and 0<z.get('pval',-1)<=1 and z.get('sebeta',-1)>0 and 0<=z.get('maf',-1)<=1
   if z and not required:reason='INVALID_OR_MISSING_REQUIRED_NATIVE_FIELD'
   row=dict(region=v['region'],rsid=v['rsid'],key=v['key'],phenocode=code,status=reason,source_sha256=qc['source_sha256'],identity_pass=identity,effect_allele=v['alt'],other_allele=v['ref'])
   for k in ['phenostring','pval','beta','sebeta','maf','maf_case','maf_control','n_case','n_control','n_sample','mlogp']:row[k]=z.get(k)
   if reason=='ESTIMATED':
    row.update(odds_ratio=math.exp(z['beta']),or_ci95_lower=math.exp(z['beta']-1.959963984540054*z['sebeta']),or_ci95_upper=math.exp(z['beta']+1.959963984540054*z['sebeta']),wald_p_diagnostic=math.erfc(abs(z['beta']/z['sebeta'])/math.sqrt(2)))
    nc,nn=z.get('n_case'),z.get('n_control');row['neff_count_formula']=4*nc*nn/(nc+nn) if finite(nc) and finite(nn) and nc>0 and nn>0 else None
   else:row.update(odds_ratio=None,or_ci95_lower=None,or_ci95_upper=None,wald_p_diagnostic=None,neff_count_formula=None)
   row['sample_size_semantics']='Native endpoint metadata/count fields; not asserted per-variant observed N';native.append(row);by[v['region'],code]=row
 slots=[]
 for s,r,a,b in SLOTS:
  x,y=by[r,a],by[r,b];ok=x['status']==y['status']=='ESTIMATED';pc=max(x['pval'],y['pval']) if ok else None
  slots.append(dict(slot_id=s,region=r,rsid=VARIANTS[r]['rsid'],sleep_code=a,outcome_code=b,sleep_status=x['status'],outcome_status=y['status'],sleep_native_p=x['pval'],outcome_native_p=y['pval'],conjunction_p=pc,bonferroni_p=min(1,4*pc) if ok else None,alpha_per_slot=ALPHA,decision=('BOTH_ASSOCIATIONS_AT_FIXED_VARIANT' if pc<=ALPHA else 'DID_NOT_MEET_FAMILY_THRESHOLD') if ok else 'NOT_ESTIMATED',classification='PHENOTYPE_SENSITIVITY',cohort_independence='UNVERIFIED',same_causal_variant='NOT_TESTED',direction_relative_discovery='NOT_A_FROZEN_ENDPOINT'))
 tsv(DEST/'fixed_variant_native_fields.tsv',native);tsv(DEST/'fixed_variant_results.tsv',slots);tsv(DEST/'fixed_variant_qc.tsv',qcs)
 print(json.dumps(slots,indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('mode',choices=['freeze','acquire','analyze']);args=p.parse_args();globals()[args.mode]()
