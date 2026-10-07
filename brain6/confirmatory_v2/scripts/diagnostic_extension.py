#!/usr/bin/env python3
"""Retrospective arithmetic and literature coordinate audit, not a new test family."""
import csv,json,math,hashlib,re
from pathlib import Path
from scipy.stats import norm
from bs4 import BeautifulSoup
from audit_evidence import ROOT,AREA,write_tsv
P=ROOT/'brain6/paper/final_package_v1';W=ROOT/'work/literature-20261007'
def rows(p):return list(csv.DictReader(p.open(),delimiter='\t'))
def number(v):
 try:x=float(v);return x if math.isfinite(x) else None
 except (ValueError,TypeError):return None

def main():
 b=[r for r in rows(P/'BRAIN6_FINAL_ALT_LOCAL_BLOCKS.tsv') if r['analysis_status']=='ESTIMATED'];replayed=[];missing=[];outside=[]
 for r in b:
  rho,var,p=map(number,(r['rho'],r['var'],r['p']))
  if rho is not None and var is not None and var>0 and p is not None:replayed.append(abs(2*norm.sf(abs(rho)/math.sqrt(var))-p))
  else:missing.append({'pair_id':r['pair_id'],'chr':r['chr'],'start':r['start'],'reason':'MISSING_OR_NONPOSITIVE_COVARIANCE_FIELD'})
  c=number(r['corr'])
  if c is not None and abs(c)>1:outside.append(r)
 d={'classification':'RETROSPECTIVE_ARITHMETIC_ONLY','estimated_blocks':len(b),'p_replayed':len(replayed),'p_missing_or_invalid':len(missing),'maximum_absolute_p_error':max(replayed),'correlation_numeric':sum(number(r['corr']) is not None for r in b),'correlation_missing':sum(number(r['corr']) is None for r in b),'correlation_outside_bounds':len(outside),'significant_outside_bounds':sum(float(r['p_fwer'])<=.05 for r in outside),'interpretation':'Covariance test and descriptive correlation are different. Out-of-range correlations cannot be interpreted as bounded population rg; no clipping or promotion.'}
 (AREA/'qc/secondary_covariance_diagnostics.json').write_text(json.dumps(d,indent=2)+'\n');write_tsv(AREA/'qc/secondary_out_of_range_correlations.tsv',outside)
 table=BeautifulSoup((W/'zu_table3.html').read_text(),'html.parser').select_one('table');pub=[];section=''
 for tr in table.select('tr'):
  v=[c.get_text(' ',strip=True) for c in tr.find_all(['td','th'],recursive=False)]
  if v and 'PLACO for' in v[0]:
   section=v[0];v=v[1:]
  if len(v)>=4 and re.fullmatch(r'rs\d+',v[1]):
   pub.append({'phenotype_section':section,'locus':v[0],'rsid':v[1],'chr':v[2],'pos':v[3],'paper_PP3':v[-3] if len(v)>7 else '', 'paper_PP4':v[-2] if len(v)>7 else '', 'source_url':'https://www.nature.com/articles/s41398-026-04166-4/tables/3','source_sha256':hashlib.sha256((W/'zu_table3.html').read_bytes()).hexdigest(),'build_status':'PAPER_BUILD_NOT_EXPLICIT_IN_RETRIEVED_HTML;rs77960_GRCh37_coordinate_independently_verified','alleles':'NOT_REPORTED_IN_THIS_TABLE'})
 write_tsv(AREA/'review/zu_table3_derived_loci.tsv',pub)
 cross=[]
 for r in rows(P/'BRAIN6_FINAL_CANDIDATES.tsv'):
  ch,st,en=map(int,re.search(r'_chr(\d+)_(\d+)_(\d+)$',r['candidate_locus_id']).groups())
  hits=[x for x in pub if 'Insomnia' in x['phenotype_section'] and r['pair_id']=='insomnia__adhd' and int(x['chr'])==ch and st<=int(x['pos'])<=en]
  cross.append({'candidate_locus_id':r['candidate_locus_id'],'pair_id':r['pair_id'],'interval_grch37':f'chr{ch}:{st}-{en}','Brain6_lead':r['lead_variants'],'prior_leads_in_interval':';'.join(x['rsid'] for x in hits),'exact_lead_match':any(x['rsid'] in r['lead_variants'].split(';') for x in hits),'interpretation':'KNOWN_SHARED_REGION_COORDINATE_CONSISTENCY;NOT_CAUSAL_VARIANT_IDENTITY' if hits else 'NO_MATCH_IN_THIS_TABLE;NOT_EVIDENCE_OF_NOVELTY','build_caveat':'Published coordinates compared provisionally; rs77960 mapping verifies chr5 example, others need complete source build attestation'})
 write_tsv(AREA/'review/locus_novelty_crosswalk.tsv',cross)
 focused={}
 for name in ['focused_sleep_gwas','focused_sleep_psychiatric']:
  dd=json.loads((W/(name+'.json')).read_text())
  for r in dd['resultList']['result']:
   key=(r.get('source'),r.get('id'));focused[key]=dict(source=r.get('source'),id=r.get('id'),title=r.get('title'),doi=r.get('doi',''),first_publication_date=r.get('firstPublicationDate',''),retrieval_query=name,screening_status='TITLE_METADATA_SCREENED',fulltext_status='NOT_REVIEWED_UNLESS_LISTED_IN_NOVELTY_AUDIT')
 write_tsv(AREA/'review/focused_search_records.tsv',list(focused.values()))
 cats=[]
 for p in sorted(W.glob('catalog_*.json')):
  dd=json.loads(p.read_text())
  if 'page' not in dd:continue
  for r in dd.get('_embedded',{}).get('associations',[]):
   cats.append({'query_object':p.name,'association_id':r['association_id'],'rs_alleles':';'.join(r.get('snp_effect_allele',[])),'GRCh38_locations':';'.join(r.get('locations',[])),'reported_trait':';'.join(r.get('reported_trait',[])),'accession':r['accession_id'],'pmid':r.get('pubmed_id',''),'p':r['p_value'],'source_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'interpretation':'Association context; not independent two-trait replication; Catalog release GRCh38.p14'})
 write_tsv(AREA/'review/catalog_association_context.tsv',cats)
 print(json.dumps(d));print('prior ADHD-insomnia region overlaps',sum(bool(x['prior_leads_in_interval']) for x in cross),'of',sum(x['pair_id']=='insomnia__adhd' for x in cross))
if __name__=='__main__':main()
