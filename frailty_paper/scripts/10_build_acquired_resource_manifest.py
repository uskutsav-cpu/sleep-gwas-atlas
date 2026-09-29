#!/usr/bin/env python3
"""Write a checksum-backed inventory of verified frailty and PubMed resources."""
from __future__ import annotations
import argparse,csv,hashlib,json,re
from datetime import datetime,timezone
from pathlib import Path
FIELDS=['resource_id','trait','resource_type','study','publication','accession','source','download_date','file','bytes','sha256','genome_build','ancestry','sample_size','cases','controls','effect_type','cohorts','ukb_overlap','finngen_overlap','license','notes']

def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''): h.update(b)
 return h.hexdigest()
def md5(p):
 h=hashlib.md5()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''): h.update(b)
 return h.hexdigest()
def load(path):
 with path.open(encoding='utf-8',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--repo',default='.');a=ap.parse_args();root=Path(a.repo).resolve();paper=root/'frailty_paper'; rows=[]
 sources={r['source_id']:r for r in load(root/'config/public_gwas_sources.tsv')}
 traits={r['trait_id']:r for r in load(root/'config/traits.tsv')}
 mapping=[('atkins_2021_frailty_index','Frailty Index','beta (continuous FI)'),('timmers_2019_parental_lifespan','Parental lifespan','beta (log-hazard protection)'),('deelen_2019_longevity_90th','Longevity, 90th percentile','beta (log odds)'),('neale_2018_left_grip_strength','Left-hand grip strength','beta (inverse-rank-normalized continuous trait)'),('bellenguez_2022_alzheimer_stage1',"Alzheimer's disease, Stage I file",'log odds ratio (beta; OR also supplied)'),('nalls_2019_parkinson_public_proxy',"Parkinson's disease (diagnosed + UK Biobank proxy)",'log odds beta')]
 for sid,label,effect in mapping:
  s=sources[sid]; tid=s['trait_ids'].split(',')[0]; t=traits[tid]; files=[x for x in s['raw_files'].split(',') if x]
  if len(files)!=1: raise ValueError(f'Expected one materialized file for {sid}: {files}')
  p=root/'data/raw'/files[0]
  if not p.is_file() or p.stat().st_size==0: raise FileNotFoundError(p)
  acc=re.search(r'GCST\d+',s['source_page_url']); accession=acc.group(0) if acc else 'NA'
  if sid=='bellenguez_2022_alzheimer_stage1': n,cases,controls='487511','85934','401577'; overlap='YES; 46,828 UK Biobank proxy cases in this Stage I file'
  elif sid=='nalls_2019_parkinson_public_proxy': n,cases,controls='482730','33674','449056'; overlap='YES; source includes 18,618 UK Biobank proxy cases'
  else: n=t['n_total'] if t['n_total'] not in ('','NA') else 'NA'; cases=t['ncase'] if t['ncase'] not in ('','NA') else 'NA'; controls=t['ncontrol'] if t['ncontrol'] not in ('','NA') else 'NA'
  if sid=='atkins_2021_frailty_index': cohorts='European-descent UK Biobank/TwinGene meta-analysis'; overlap='YES; stated in registered trait metadata'
  elif sid=='neale_2018_left_grip_strength': cohorts='UK Biobank, field 46, both sexes'; overlap='YES; UK Biobank source'
  elif sid=='bellenguez_2022_alzheimer_stage1': cohorts='European Stage I; clinically diagnosed AD and UK Biobank proxy cases';
  elif sid=='nalls_2019_parkinson_public_proxy': cohorts='European discovery with diagnosed and UK Biobank proxy cases; see source registry and publication';
  else: cohorts='Cohort list not fully resolved in local locked source registry; see source publication'; overlap='UNKNOWN; verify against source publication before overlap-sensitive analysis'
  build='GRCh37' if sid in {'atkins_2021_frailty_index','timmers_2019_parental_lifespan','deelen_2019_longevity_90th','neale_2018_left_grip_strength','nalls_2019_parkinson_public_proxy'} else 'GRCh38'
  notes=f"Local source registry: {s['source_id']}; {s['notes']}"
  rows.append({'resource_id':sid,'trait':label,'resource_type':'GWAS summary statistics','study':sid,'publication':f"PMID:{s['pmid']}" if s['pmid'] else 'Unpublished public release; no PMID in registry','accession':accession,'source':s['download_url'],'download_date':'2026-09-22 (collection log date)','file':str(p.relative_to(root)),'bytes':p.stat().st_size,'sha256':sha(p),'genome_build':build,'ancestry':s['ancestry_reported'] or 'UNKNOWN','sample_size':n,'cases':cases,'controls':controls,'effect_type':effect,'cohorts':cohorts,'ukb_overlap':overlap,'finngen_overlap':'UNKNOWN; no explicit cohort overlap statement in the local source registry','license':'Public access; reuse license not separately recorded in local source registry','notes':notes})
 # Register the frozen construct-replication endpoint only after bytes exist.
 # Source metadata remains explicit about unresolved build/effect conventions;
 # acquisition must not implicitly promote the file into harmonized analyses.
 ffs_path=paper/'data/gwas/fried_frailty_score/Frailty_GWAS_maf0.01_info0.9_hwe1e-6_FUMA.tsv'
 if ffs_path.is_file() and ffs_path.stat().st_size>0:
  rows.append({'resource_id':'ye_2023_fried_frailty_score','trait':'Fried Frailty Score (ordinal 0–5)','resource_type':'GWAS summary statistics','study':'Ye et al. 2023 UK Biobank Fried Frailty Score GWAS','publication':'PMID:36928559; DOI:10.1007/s11357-023-00771-z','accession':'GCST90295968','source':'https://figshare.com/s/6683396c68807fe4e729','download_date':'2026-09-26','file':str(ffs_path.relative_to(root)),'bytes':ffs_path.stat().st_size,'sha256':sha(ffs_path),'genome_build':'UNKNOWN; verify against this exact full-statistics release before harmonization','ancestry':'European ancestry','sample_size':'386565','cases':'NA','controls':'NA','effect_type':'File columns BETA/SE/P/OR; article reports BOLT-LMM single-variant testing for ordinal FFS. OR equals exp(BETA) in the acquired file; treat BETA as the reported effect and do not interpret OR as a separate odds-ratio estimate without source confirmation.','cohorts':'UK Biobank discovery','ukb_overlap':'YES; UK Biobank cohort; exact participant overlap with each sleep source UNKNOWN','finngen_overlap':'UNKNOWN; no explicit statement of zero overlap','license':'CC BY 4.0 (Figshare share page)','notes':'Full file from the article-linked Figshare share page; ndownloader range probe established 639246335 total bytes, and the downloaded MD5 matched the S3 ETag. Stream validation: 8,883,488 data rows plus header; exact 9-column schema; no row-width, numeric, position, p-value, SE, OR, or allele-domain anomalies. These checks establish acquisition integrity, not genome-build validation, effect-allele semantics, harmonization, cohort independence, or analysis eligibility.'})
 # Add each phenotype in the locked 12-trait sleep panel as a distinct row.
 # A shared source (accelerometer) therefore produces three trait rows while
 # still pointing to its exact registered raw file and common archive.
 panel_path=root/'config/analysis_panel.tsv'
 if panel_path.is_file():
  panel=[r for r in load(panel_path) if r.get('domain')=='sleep']
  if len(panel)!=12 or len({r['trait_id'] for r in panel})!=12:
   raise ValueError(f'Expected 12 unique sleep traits in locked atlas panel, found {len(panel)}')
  for trait in panel:
   sid=trait['source_id']; s=sources.get(sid)
   if not s or s.get('access')!='PUBLIC' or s.get('acquisition_status')!='SOURCE_VERIFIED_PUBLIC':
    raise ValueError(f'Sleep panel source is not registered as verified public: {sid}')
   raw=root/'data/raw'/trait['raw_file']
   archive=root/'data/raw/.archives'/s['archive_name']
   # Skip genuinely unacquired inputs so the manifest remains usable during
   # resumable collection. Never register a partial archive or missing raw file.
   if not archive.is_file() or not raw.is_file() or archive.stat().st_size!=int(s['archive_bytes']):
    continue
   if sha(archive)!=s['archive_sha256'].lower():
    raise ValueError(f'Registered archive checksum/size mismatch for {sid}: {archive}')
   if raw.stat().st_size==0:
    raise ValueError(f'Empty materialized sleep GWAS file: {raw}')
   build={'hg19':'GRCh37','hg38':'GRCh38'}.get(trait.get('build','').lower(),'UNKNOWN')
   note=trait.get('source_note','')
   ukb='YES; source note identifies UK Biobank' if 'uk biobank' in note.lower() or 'ukb' in note.lower() else 'UNKNOWN; source-specific cohort overlap unresolved'
   finngen='YES; FinnGen source' if 'finngen' in sid.lower() else 'UNKNOWN; source-specific cohort overlap unresolved'
   rows.append({'resource_id':'sleep_panel_'+trait['trait_id'],'trait':trait['label'],'resource_type':'GWAS summary statistics','study':trait.get('dataset_version') or sid,'publication':f"PMID:{trait['pmid']}" if trait.get('pmid') else 'Publication identifier not recorded in locked panel','accession':'NA; no GWAS Catalog accession recorded in locked panel','source':s['download_url'],'download_date':'2026-09-23 (registered public archive collection)','file':str(raw.relative_to(root)),'bytes':raw.stat().st_size,'sha256':sha(raw),'genome_build':build,'ancestry':trait.get('ancestry') or s.get('ancestry_reported') or 'UNKNOWN','sample_size':trait.get('n_total') or 'UNKNOWN','cases':trait.get('ncase') or 'NA','controls':trait.get('ncontrol') or 'NA','effect_type':f"Beta/effect scale as in source; phenotype type {trait.get('type','UNKNOWN')}",'cohorts':note or 'Not fully specified in locked panel; consult source study','ukb_overlap':ukb,'finngen_overlap':finngen,'license':'Public source; reuse license not separately recorded in source registry','notes':f"Locked sleep panel trait {trait['trait_id']}; source_id={sid}; registered archive size/SHA-256 verified; panel build={trait.get('build','UNKNOWN')}; source provenance={note}"})
 pub=paper/'review/pubmed_source_files.tsv'; pubrows=load(pub)
 qlabel={'primary_sleep_frailty':'Primary sleep–frailty PubMed search','genetic_sleep_frailty':'Genetic sleep–frailty PubMed search','frailty_neighborhood':'Secondary frailty/aging-neighborhood PubMed search'}
 for x in pubrows:
  rel='frailty_paper/review/'+x['source_file']
  source_path=root/rel
  downloaded=datetime.fromtimestamp(source_path.stat().st_mtime,timezone.utc).date().isoformat()
  q=x['source_file'].split('/')[1] if '/' in x['source_file'] else 'UNKNOWN'
  rows.append({'resource_id':'pubmed_xml_'+hashlib.sha256(rel.encode()).hexdigest()[:16],'trait':qlabel.get(q,q),'resource_type':'PubMed XML search result batch','study':q,'publication':'NA; bibliographic search records','accession':'NA','source':'NCBI Entrez E-Utilities','download_date':f'{downloaded} (XML file timestamp, UTC; see search/reconciliation logs)','file':rel,'bytes':x['bytes'],'sha256':x['sha256'],'genome_build':'NA','ancestry':'NA','sample_size':'NA','cases':'NA','controls':'NA','effect_type':'NA','cohorts':'NA','ukb_overlap':'NA','finngen_overlap':'NA','license':'Public bibliographic metadata; record-level reuse terms not adjudicated','notes':f"{x['record_count']} PubMed records; exact query/count history in frailty_paper/review/pubmed/search_log.tsv and preserved dated logs; supplemental index updates in frailty_paper/review/pubmed_update_log.tsv; PMID/DOI/title dedup audit in review/all_records.tsv and deduplicated_records.tsv."})
 indexes=[
  ('eqtl_catalogue_dataset_metadata_r7','eQTL Catalogue dataset metadata index','eQTL Catalogue dataset metadata table',paper/'data/qtl/eqtl_catalogue_dataset_metadata_r7.tsv','https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/data_tables/dataset_metadata_r7.tsv','Public repository metadata; downstream dataset terms not adjudicated','758 metadata rows; index only, no eQTL summary statistics acquired.'),
  ('brainscope_key_resource_links','BrainSCOPE key-resource link index','BrainSCOPE key-resource page',paper/'data/single_cell/brainscope_key_resource_links.json','https://brainscope.gersteinlab.org/key_resource_files.html','Public resource index; underlying resources have individual terms','29 links captured from the page; no underlying expression, QTL, or sample files acquired.'),
  ('brainscope_sample_metadata','BrainSCOPE DLPFC sample metadata','BrainSCOPE PEC2 sample metadata',paper/'data/single_cell/brainscope_key_files/PEC2_sample_metadata.txt','https://brainscope.gersteinlab.org/data/sample_metadata/PEC2_sample_metadata.txt','Public resource metadata; source-specific terms should be checked before downstream use','Public sample metadata only; no expression matrices or QTL statistics acquired.'),
  ('brainscope_eqtl_links','BrainSCOPE cell-type eQTL link index','BrainSCOPE cell-type eQTL page',paper/'data/single_cell/brainscope_eqtl_links.json','https://brainscope.gersteinlab.org/output-sig-eQTL.html','Public resource index; underlying resources have individual terms','23 links captured from the page; no QTL archive or cell-type summary statistics acquired.')]
 for rid,label,study,p,source,license,notes in indexes:
  if not p.is_file() or p.stat().st_size==0: continue
  rows.append({'resource_id':rid,'trait':label,'resource_type':'Public resource index / metadata','study':study,'publication':'NA; official resource page','accession':'NA','source':source,'download_date':'2026-09-22 (source-index collection log date)','file':str(p.relative_to(root)),'bytes':p.stat().st_size,'sha256':sha(p),'genome_build':'NA','ancestry':'NA','sample_size':'NA','cases':'NA','controls':'NA','effect_type':'NA','cohorts':'NA','ukb_overlap':'NA','finngen_overlap':'NA','license':license,'notes':notes})
 component_dir=paper/'data/gwas/physical_frailty_components'
 metadata_path=component_dir/'zenodo_metadata.json'
 component_counts={
  'GWAS_Weight_loss_maf0.01_info0.9.txt.gz':('Weight loss','61378','340745'),
  'GWAS_Exhaustion_maf0.01_info0.9.txt.gz':('Exhaustion','47139','349680'),
  'GWAS_Low_Physical_activity_maf0.01_info0.9.txt.gz':('Low physical activity','34406','371518'),
  'GWAS_Slow_Walking_speed_maf0.01_info0.9.txt.gz':('Slow walking speed','30846','375531'),
  'GWAS_Low_Grip_strength_maf0.01_info0.9.txt.gz':('Low grip strength','56262','350718'),
 }
 if metadata_path.is_file():
  meta=json.loads(metadata_path.read_text(encoding='utf-8'))
  meta_files={x.get('key'):x for x in meta.get('files',[])}
 for name,(trait,cases,controls) in component_counts.items():
   source_meta=meta_files.get(name); p=component_dir/name
   if not source_meta or not p.is_file() or p.stat().st_size!=int(source_meta['size']): continue
   checksum=source_meta.get('checksum','')
   if not checksum.startswith('md5:') or md5(p)!=checksum.split(':',1)[1].lower(): continue
   doi=meta.get('metadata',{}).get('doi','10.5281/zenodo.14011550')
   source=(source_meta.get('links') or {}).get('self','https://zenodo.org/api/records/14011550/files/'+name+'/content')
   rows.append({'resource_id':'zenodo_14011550_'+re.sub(r'[^a-z0-9]+','_',trait.lower()).strip('_'),'trait':trait+' physical frailty component','resource_type':'GWAS summary statistics','study':'Zenodo record 14011550','publication':'Chen et al.; deposit cites Physical Frailty, Inflammation, Polygenic Risk Score, and Primary Prevention of Heart Failure in Patients with Cardiovascular Diseases','accession':doi,'source':source,'download_date':'2026-09-22 (local file timestamp)','file':str(p.relative_to(root)),'bytes':p.stat().st_size,'sha256':sha(p),'genome_build':'UNKNOWN; not stated in deposit README','ancestry':'UNKNOWN; not stated in deposit README','sample_size':str(int(cases)+int(controls)),'cases':cases,'controls':controls,'effect_type':'BETA supplied; scale/model to be checked against study methods','cohorts':'UK Biobank (stated in Zenodo record description)','ukb_overlap':'YES; UK Biobank cohort stated in Zenodo record','finngen_overlap':'UNKNOWN','license':meta.get('metadata',{}).get('license',{}).get('id','UNKNOWN'),'notes':f"File size and Zenodo MD5 verified. Counts from {component_dir.name}/0_Readme.txt; methods, ancestry and build still require primary-study verification before harmonization."})
 health_dir=paper/'data/gwas/healthspan'
 health_meta_path=health_dir/'zenodo_metadata.json'
 if health_meta_path.is_file():
  hmeta=json.loads(health_meta_path.read_text(encoding='utf-8'))
  item=next((x for x in hmeta.get('files',[]) if x.get('key')=='healthspan_summary.csv.gz'),None)
  hp=health_dir/'healthspan_summary.csv.gz'
  if item and hp.is_file() and hp.stat().st_size==int(item.get('size',-1)):
   checksum=item.get('checksum','')
   if checksum.startswith('md5:') and md5(hp)==checksum.split(':',1)[1].lower():
    source=(item.get('links') or {}).get('content') or (item.get('links') or {}).get('self','https://zenodo.org/api/records/1302861/files/healthspan_summary.csv.gz/content')
    rows.append({'resource_id':'zenodo_1302861_healthspan','trait':'Human healthspan','resource_type':'GWAS summary statistics','study':'Zenin et al. 2019','publication':'Zenin et al. 2019, Communications Biology, doi:10.1038/s42003-019-0290-0','accession':hmeta.get('metadata',{}).get('doi','10.5281/zenodo.1302861'),'source':source,'download_date':'2026-09-22 (local file timestamp)','file':str(hp.relative_to(root)),'bytes':hp.stat().st_size,'sha256':sha(hp),'genome_build':'GRCh37 (Zenodo record description)','ancestry':'Genetically Caucasian, British; as described in Zenodo record','sample_size':'300447','cases':'NA','controls':'NA','effect_type':'Beta, SE, Z and -log10(P) supplied; verify header and trait coding during QC','cohorts':'UK Biobank','ukb_overlap':'YES; UK Biobank','finngen_overlap':'NO stated; verify against analysis cohort manifest','license':hmeta.get('metadata',{}).get('license',{}).get('id','UNKNOWN'),'notes':'File size and Zenodo MD5 verified. Sample size, ancestry and build from the official Zenodo record description; source publication and data-use conditions recorded in record metadata.'})
 latent_index=paper/'manifests/latent_frailty_accession_index.tsv'
 if latent_index.is_file():
  for x in load(latent_index):
   if x.get('local_acquisition_status')!='ACQUIRED_MD5_VERIFIED': continue
   acc=x['accession']; lp=paper/'data/gwas/latent_frailty_catalog'/acc/(acc+'.tsv')
   checks=paper/'data/gwas/latent_frailty_catalog'/acc/'md5sum.txt'
   expected=''
   if checks.is_file():
    match=re.search(rf'^([0-9a-fA-F]{{32}})\s+\*?{re.escape(acc+".tsv")}\s*$',checks.read_text(encoding='ascii'),re.M)
    expected=match.group(1).lower() if match else ''
   if not lp.is_file() or lp.stat().st_size==0 or not expected or md5(lp)!=expected: continue
   label=x.get('phenotype_label','').strip()
   if x.get('label_status')!='VERIFIED' or not label:
    label=f"Unresolved phenotype, GWAS Catalog accession {acc}"
   rows.append({'resource_id':'gwas_catalog_'+acc,'trait':label,'resource_type':'GWAS summary statistics','study':'Foote et al. multivariate frailty GWAS; GWAS Catalog '+acc,'publication':'Foote et al. 2025, Nature Genetics, doi:10.1038/s41588-025-02269-0','accession':acc,'source':x.get('ftp_tsv_url',''),'download_date':'2026-09-22 (local file timestamp)','file':str(lp.relative_to(root)),'bytes':lp.stat().st_size,'sha256':sha(lp),'genome_build':x.get('genome_assembly','UNKNOWN'),'ancestry':x.get('sample_ancestry','UNKNOWN'),'sample_size':x.get('sample_size','UNKNOWN'),'cases':'NA','controls':'NA','effect_type':'Beta; consult accession metadata for N definition and Q metric','cohorts':'Not resolved by accession metadata alone','ukb_overlap':'UNKNOWN; cohort composition must be confirmed from study methods','finngen_overlap':'UNKNOWN','license':x.get('license','UNKNOWN'),'notes':f"Official Catalog md5sum.txt verified; accession-to-phenotype mapping status: {x.get('label_status','UNRESOLVED')}. Retain accession-only identity when unresolved; do not infer a factor label."})
 # Register retained HLMA article supplements and portal metadata snapshots.
 # These are metadata inputs only; no processed muscle archives are registered
 # as acquired from the listing or its link targets.
 existing_ids={row.get('resource_id') for row in rows}
 hlma_metadata=paper/'manifests/HLMA_PMC_Article_Datasets_PMC11062927.1.json'
 supplement_specs=[
  ('hlma_supplementary_table_1_donor_information','Donor-level atlas metadata','MOESM3_ESM.xlsx','HLMA_Supplementary_Table_1_Donor_Information.xlsx','Supplementary Table 1','European and Asian Chinese cohorts','31 donor rows','European (18); Asian Chinese (13)','PMC Article Datasets media entry MD5 verified; XLSX package valid; one sheet, 31 donor rows, 12 columns. Donor-level demographic/clinical descriptions only; do not infer identity.'),
  ('hlma_supplementary_table_2_library_information','Atlas library and assay metadata','MOESM4_ESM.xlsx','HLMA_Supplementary_Table_2_Library_Information.xlsx','Supplementary Table 2','NA','190 library records across 3 assays','scRNA-seq 38; snRNA-seq 106; snATAC-seq 46','PMC Article Datasets media entry MD5 verified; XLSX package valid; one sheet with modality-specific headers; cell/nucleus-count sums: scRNA-seq 79,649; snRNA-seq 212,774; snATAC-seq 95,021. scRNA donor code YM5 is absent from Supplementary Table 1; crosswalk unresolved.'),
 ]
 if hlma_metadata.is_file():
  meta=json.loads(hlma_metadata.read_text(encoding='utf-8'))
  media=meta.get('media_urls',[])
  for rid,trait,suffix,filename,table,ancestry,n,cohorts,notes in supplement_specs:
   if rid in existing_ids: continue
   p=paper/'data/single_cell/hlma_supplementary'/filename
   url=next((u for u in media if suffix in u),None)
   match=re.search(r'[?&]md5=([0-9a-fA-F]{32})',url or '')
   if not p.is_file() or not url or not match or md5(p)!=match.group(1).lower():
    continue
   rows.append({'resource_id':rid,'trait':trait,'resource_type':'Supplementary donor metadata' if 'table_1' in rid else 'Supplementary library metadata','study':'HLMA_2024','publication':f"PMID:{meta['pmid']}; DOI:{meta['doi']}",'accession':f"{meta['pmcid']}.{meta['version']}; {table}",'source':url,'download_date':'2026-09-23','file':str(p.relative_to(root)),'bytes':p.stat().st_size,'sha256':sha(p),'genome_build':'NA','ancestry':ancestry,'sample_size':n,'cases':'NA','controls':'NA','effect_type':'NA','cohorts':cohorts,'ukb_overlap':'NA','finngen_overlap':'NA','license':meta.get('license_code','UNKNOWN')+' (PMC article-version metadata)','notes':notes})
 # Archive the publisher's public HFRS supplementary hit tables when present.
 # These are thresholded reported-variant tables, not genome-wide summary stats.
 hfrs_supplement=paper/'data/gwas/hfrs_mak_2025/mak_2025_hfrs_supplementary_tables.xlsx'
 if hfrs_supplement.is_file():
  expected_sha='123c340cf063e2ef618b8ef543dedea6b7c8ef69008488f01027eedbffc04625'
  observed_sha=sha(hfrs_supplement)
  if observed_sha!=expected_sha:
   raise ValueError(f'Unexpected publisher HFRS supplement SHA-256: {observed_sha}')
  rows.append({'resource_id':'mak_2025_hfrs_published_supplementary_hit_tables','trait':'Hospital Frailty Risk Score (HFRS), with- and without-dementia hit tables','resource_type':'Publisher supplementary GWAS hit tables (not full summary statistics)','study':'Mak et al. 2025 FinnGen R12 HFRS GWAS','publication':'Mak et al. 2025, Nature Aging 5:1589-1600, doi:10.1038/s43587-025-00925-y','accession':'Nature Aging supplementary information; Supplementary Tables 1-19','source':'https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs43587-025-00925-y/MediaObjects/43587_2025_925_MOESM1_ESM.xlsx','download_date':'2026-09-24','file':str(hfrs_supplement.relative_to(root)),'bytes':hfrs_supplement.stat().st_size,'sha256':observed_sha,'genome_build':'Not stated in workbook; confirm against exact source release before locus harmonization','ancestry':'FinnGen Finnish discovery; UK Biobank white British replication as described in article','sample_size':'FinnGen discovery N=500737; UK Biobank replication N=407463','cases':'NA','controls':'NA','effect_type':'FinnGen beta/SE/P and UKB beta/SE/P plus METAL Z/P for reported FinnGen P<5e-8 hits only','cohorts':'FinnGen R12 discovery; UK Biobank replication; combined METAL meta-analysis','ukb_overlap':'YES; UK Biobank replication cohort; exact overlap with project sleep/frailty inputs unresolved','finngen_overlap':'YES; FinnGen discovery cohort; exact overlap with project sources unresolved','license':'CC-BY-4.0 (article rights statement; cite the article)','notes':'Official publisher workbook, package validated. ST1 has 1,588 nonempty HFRS result rows and ST2 has 492 nonempty HFRS-without-dementia result rows; both contain reported P<5e-8 hits and UKB/METAL columns, not full genome-wide summary statistics. Do not use for LDSC/LAVA or claim independent UKB replication. See frailty_paper/analysis/hfrs_published_supplement_audit_2026-09-24.md.'})
 portal_snapshots=sorted((paper/'manifests').glob('hlma_portal_download_metadata_*.json'))
 for p in portal_snapshots:
  snapshot=json.loads(p.read_text(encoding='utf-8'))
  rid='hlma_portal_download_metadata_'+p.stem.rsplit('_',1)[-1]
  if rid in existing_ids: continue
  portal_rows=sum(len(section.get('data',[])) for section in snapshot.get('data',[]))
  if not snapshot.get('data') or portal_rows==0:
   raise ValueError(f'Empty HLMA portal listing snapshot: {p}')
  rows.append({'resource_id':rid,'trait':'HLMA processed-object portal listing metadata','resource_type':'Public resource index / metadata','study':'HLMA_2024 official download listing','publication':'PMID:38649488; DOI:10.1038/s41586-024-07348-6','accession':'PRJCA017356; OMIX IDs are recorded in the listing links','source':'https://db.genomics.cn/cdcp_hlma/api/download/get_download','download_date':p.stem.rsplit('_',1)[-1],'file':str(p.relative_to(root)),'bytes':p.stat().st_size,'sha256':sha(p),'genome_build':'NA','ancestry':'NA','sample_size':f'{portal_rows} portal listing rows; metadata only','cases':'NA','controls':'NA','effect_type':'NA','cohorts':'HLMA atlas; file-level donor linkage unresolved for some rows','ukb_overlap':'NA','finngen_overlap':'NA','license':'Portal listing metadata only; no object-specific reuse license inferred','notes':f'{portal_rows} listed object/fragment rows with public portal file_path metadata. Link identifiers do not establish archive byte identity, archive-member identity, donor crosswalk, coordinate build, or reuse terms; see frailty_paper/analysis/hlma_portal_omix_crosswalk_2026-09-24.md.'})
 out=paper/'manifests/all_acquired_resources.tsv';out.parent.mkdir(parents=True,exist_ok=True)
 with out.open('w',newline='',encoding='utf-8') as f:
  w=csv.DictWriter(f,fieldnames=FIELDS,delimiter='\t',lineterminator='\r\n');w.writeheader();w.writerows(rows)
 gwas_count=sum(1 for row in rows if row['resource_type']=='GWAS summary statistics')
 print(f'Wrote {out}: {len(rows)} rows ({gwas_count} GWAS files; {len(pubrows)} PubMed XML batches)')
if __name__=='__main__':main()
