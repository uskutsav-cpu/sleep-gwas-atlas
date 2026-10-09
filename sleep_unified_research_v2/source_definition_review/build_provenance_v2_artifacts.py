"""Serialize result-free public-source review evidence; performs no network/GWAS work."""
from pathlib import Path
from datetime import datetime, timezone
import csv
import hashlib
import json
import re

REVIEW = Path(__file__).resolve().parent
CACHE = Path('/Volumes/Extreme SSD/sleep-unified-research-v1/recovery-2026-10-09/source_definition_review_v2')
INITIAL = CACHE / 'initial_cache'

CATALOG = 'https://www.ebi.ac.uk/gwas/rest/api/v2/studies/GCST90475826'
YAML = 'https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90475001-GCST90476000/GCST90475826/GCST90475826.tsv.gz-meta.yaml'
GIA = 'https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/analysis.cgi?study_id=phs002453.v1.p1&pha=10373'
DICT = 'https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/document.cgi?phd=8759&study_id=phs002453.v1.p1'
README = 'https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/GetPdf.cgi?id=phd008760.1'
SUPP = 'https://pmc-oa-opendata.s3.amazonaws.com/PMC12857194.1/NIHMS2091849-supplement-Supplemental_Text.pdf'
SUPPMETA = 'https://pmc-oa-opendata.s3.amazonaws.com/PMC12857194.1/PMC12857194.1.json'
CIPHER = 'https://phenomics.va.ornl.gov/web/api/phenotype/14926'
CIPHER_UI = 'https://phenomics.va.ornl.gov/web/cipher/phenotype-viewer/details?id=14926'
AUTHOR = 'https://zenodo.org/records/15595560/files/ady_response_huffman.v2.pdf?download=1'
PHEWAS = 'https://github.com/PheWAS/PheWAS/tree/55dd1c24e228851922400cfba8d7db474565ccc7'

def dump(name, data):
    (REVIEW/name).write_text(json.dumps(data, indent=2, ensure_ascii=False)+'\n')

def digest(p, algorithm='sha256'):
    return getattr(hashlib, algorithm)(p.read_bytes()).hexdigest()

def tsv(name, fields, rows):
    with (REVIEW/name).open('w', newline='') as f:
        w = csv.DictWriter(f, fields, delimiter='\t')
        w.writeheader()
        w.writerows(rows)

j = json.loads((CACHE/'cipher_insomnia_metadata.json').read_text())
catalog = json.loads((INITIAL/'catalog.json').read_text())
pmc = json.loads((CACHE/'pmc_article_metadata.json').read_text())
assert j['id'] == 14926 and j['fullName'] == 'Insomnia (gwPheWAS)'
assert 'Phe_327_4' in [r['keyword'] for r in j['keywords']]
assert catalog['accession_id'] == 'GCST90475826'
assert catalog['initial_sample_size'] == '78,566 European ancestry cases, 329,572 European ancestry controls'
assert pmc['pmid'] == 39024449 and pmc['doi'] == '10.1126/science.adj1182'
assert (CACHE/'gia_supplement.pdf').read_bytes().startswith(b'%PDF-')
assert digest(CACHE/'gia_supplement.pdf', 'md5') == 'afe2d8e36ab1a5213f3ec32a14c19a8a'

code_lists = {}
mapping = []
for group in j['algorithm']['assocCodes']:
    vocabulary = {460: 'ICD9CM', 461: 'ICD10CM'}[group['codeType']]
    code_lists[vocabulary] = [r['code'] for r in group['codes']]
    for r in group['codes']:
        mapping.append({'accession':'GCST90475826','phecode':'327.4','cipher_phenotype_id':14926,
                        'vocabulary':vocabulary,'code':r['code'],'evidence_status':'VERIFIED_OFFICIAL_CIPHER_LIST',
                        'source_url':CIPHER,'cipher_revision':j['revision'],
                        'mapping_release_caveat':'current official gwPheWAS record; no executable historical extraction artifact'})
assert len(mapping) == 18 and len(code_lists['ICD9CM']) == 8 and len(code_lists['ICD10CM']) == 10
tsv('provenance_v2_mvp_insomnia_icd_mapping.tsv', list(mapping[0]), mapping)

claims = [
 ('exact_catalog_identity','VERIFIED_ACCESSION_METADATA','MVP; insomnia PheCode327.4; EUR78566cases/329572controls; N408138; PMID39024449',CATALOG,'Catalog metadata','No body rows or outcomes read'),
 ('exact_gia_identity','VERIFIED_ACCESSION_METADATA','pha010373.1; Phe_327_4.EUR.GIA; exact same EUR cases/controls/N',GIA,'Analysis name/description','Links release identity, not implementation source code'),
 ('gia_hare_separation','VERIFIED_DICTIONARY_METADATA','GIA EUR78566/329572; separate HARE EUR80037/334480',DICT,'GIA and HARE Phe_327_4 rows','Dictionary cached body hash; redundant body discarded during disk shortage'),
 ('hare_only_document','VERIFIED_EXPLICIT_HARE_SCOPE','README detailed appendix explicitly HARE; EHR September2019',README,'PDF pages1–2','Cannot alone certify GIA'),
 ('gia_publication_ancestry','VERIFIED_PRIMARY_FINAL_GIA','1000G reference PCs; random forest assignment probability >50%; AFR/AMR/EAS/EUR analysed',SUPP,'PDF page2','SAS membership described but insufficient cases for these analyses'),
 ('ehr_end_horizon','VERIFIED_PRIMARY_FINAL_GIA','EHR follow-up through September2020',SUPP,'PDF page2','No exact start date/minimum observation interval'),
 ('phenotype_resource_link','VERIFIED_PRIMARY_FINAL_GIA','Initial phenotypes from CIPHER; ICD clinical phenotypes mapped to1876PheCodes',SUPP,'PDF page3','Links resource methods to exact phenotype record through resource/trait/count chain'),
 ('exact_cipher_identity','VERIFIED_OFFICIAL_PHENOTYPE_METADATA','ID14926 Insomnia(gwPheWAS); keywordPhe_327_4; MVP enrollees; revision3',CIPHER,'fullName/description/keywords','CIPHER record does not itself contain GCST or exact EUR counts'),
 ('icd_code_list','VERIFIED_OFFICIAL_CIPHER_LIST','8 ICD9CM codes plus10 ICD10CM codes; serialized separately',CIPHER,'algorithm.assocCodes','Do not substitute a default PheWAS map; MVP reports manual mapping'),
 ('case_occurrences','VERIFIED_PRIMARY_GIA_AND_CIPHER','At least2 target-mapped ICD instances',SUPP+' ; '+CIPHER,'PDF page3; algorithm.algorithmDesc','Distinct dates/visits and minimum spacing not specified'),
 ('control_occurrences','VERIFIED_PRIMARY_GIA_AND_CIPHER','Zero target-mapped ICD instances',SUPP+' ; '+CIPHER,'PDF page3; algorithm.algorithmDesc','Do not infer healthy controls or absence of other sleep disorders'),
 ('single_occurrence','INFERRED_FROM_EXPLICIT_RULE','One instance satisfies neither >=2 case nor0 control criterion',SUPP+' ; '+CIPHER,'Logical consequence of rule','No executable phenotyping implementation audited'),
 ('related_phecode_exclusions','VERIFIED_EXPLICIT_NONE','PheCode exclusions were not applied',CIPHER,'algorithm.algorithmDesc','Standard PheWAS default exclusions must not be imported'),
 ('manual_icd_mapping','VERIFIED_OFFICIAL_CIPHER_STATEMENT','MVP performed its own manual ICD9-to-ICD10 mapping',CIPHER,'algorithm.algorithmDesc','versionInfo1.2 is not certified upstream mapping/software version'),
 ('distinct_dates','UNRESOLVED_NOT_STATED','Instance counting versus distinct-date counting unresolved',CIPHER+' ; '+SUPP,'No explicit rule or attached executable code','Do not claim >=2 separate dates or visits'),
 ('lookback_start','UNRESOLVED_NOT_STATED','No minimum lookback or encounter density; CIPHER dataUsedStart/dataUsedEnd=null',CIPHER,'algorithm date fields','Publication supplies September2020 horizon only'),
 ('minimum_code_interval','UNRESOLVED_NOT_STATED','No minimum interval or chronicity requirement specified',CIPHER+' ; '+SUPP,'No specified rule','No invented chronic insomnia criterion'),
 ('historical_implementation','UNRESOLVED_IMPLEMENTATION','Current CIPHER list is documented; release-pinned executable historical mapping/extraction absent',CIPHER,'revision3; codeSamples empty','Source chain is primary resource/trait/count evidence, not individual-level verification'),
 ('validation','NO_VALIDATION_REPORTED','CIPHER algorithm validated=false; validationDescription=null',CIPHER,'algorithm validation fields','Not evidence that a diagnosis itself is clinically invalid'),
 ('reference_map','REFERENCE_ONLY_NOT_MVP_SUBSTITUTE','Maintainer map+rollup give matching18 insomnia codes; configurable aggregation/default exclusions',PHEWAS,'Pinned code/data/docs','MVP manual mapping and explicit no-exclusion rule remain authoritative'),
 ('effect_allele','VERIFIED_AUTHOR_EXPLANATION_NOT_INDEPENDENT_AUDIT','OR/effect tied to explicitEA; REF/ALT tested alleles may differ from genomic-reference roles',AUTHOR,'PDF pages1–2','Does not establish EBI conversion or SE scale'),
 ('pheweb_orientation','VERIFIED_AUTHOR_EXPLANATION_NOT_INDEPENDENT_AUDIT','Some PheWeb views use minor-allele orientation and may differ from dbGaPEA',AUTHOR,'PDF page2','Do not mix plot orientation with raw-source effect allele'),
 ('mac_filter_discrepancy','UNRESOLVED_RELEASE_OR_STAGE_DISCREPANCY','Exact dbGaP MAC<30 filter; supplement genotypeQC MAC<20 removal and GWAS MAC>40 inclusion',GIA+' ; '+SUPP,'Exact analysis methods; PDF pages2,4','Could differ by stage/release; public evidence does not resolve explanation'),
 ('rights_separation','VERIFIED_DISTINCT_METADATA_LICENSES','GWAS Catalog accession EBIterms; PMC source-article JSON CC BY; no GWAS-license transfer implied',CATALOG+' ; '+SUPPMETA,'terms_of_license; license_code','Article/supplement rights are not summary-statistic rights'),
]
rows=[dict(claim_id=c,status=s,finding=f,source_urls=u,source_locator=l,limitation=q) for c,s,f,u,l,q in claims]
tsv('provenance_v2_definition_evidence.tsv',list(rows[0]),rows)

contract = {
 'review_version':'v2','gwas_accession':'GCST90475826','dbgap_analysis':'pha010373.1',
 'analysis_name':'Phe_327_4.EUR.GIA','publication_pmid':39024449,
 'cohort':'MVP','ancestry':'European GIA','n_cases':78566,'n_controls':329572,'n_total':408138,
 'ascertainment':'clinical EHR ICD-coded insomnia; related phenotype to core UKB frequent self-report insomnia',
 'code_lists':code_lists,'code_list_authority':CIPHER,
 'case_rule':{'mapped_icd_instances_minimum':2,'status':'VERIFIED_PRIMARY_GIA_AND_CIPHER'},
 'control_rule':{'mapped_icd_instances':0,'related_phecode_exclusions_applied':False,'status':'VERIFIED_PRIMARY_GIA_AND_CIPHER'},
 'single_instance':{'eligible_as_stated_case_or_control':False,'status':'INFERRED_FROM_EXPLICIT_RULE'},
 'ehr_end_horizon':{'value':'September 2020','status':'VERIFIED_FINAL_PUBLICATION','source':SUPP,'pdf_page':2},
 'unresolved': ['distinct-date versus repeated-code instance counting','minimum interval between code instances',
                'start lookback date and minimum observation duration','historical executable mapping/extraction artifact',
                'accession MAC<30 versus final-publication GWAS MAC>40 discrepancy','effect SE scale and EBI conversion independent validation'],
 'cipher_snapshot':{'id':14926,'revision':j['revision'],'versionInfo':j['versionInfo'],
                    'created':j['created'],'lastModified':j['lastModified'],'algorithmCreated':j['algorithm']['algorithmCreated'],
                    'revisionComment':j['revisionComment'],'majorRevision':j['majorRevision'],
                    'validated':j['algorithm']['validated'],'dataUsedStart':None,'dataUsedEnd':None},
 'source_linkage_type':'PRIMARY_RESOURCE_TRAIT_COUNT_CHAIN_NOT_PARTICIPANT_IMPLEMENTATION_AUDIT',
 'source_linkage': [CATALOG,GIA,DICT,SUPP,CIPHER],
 'effect_allele_contract':{'rule':'Use explicit effect allele and matching tested allele pair; do not equate raw ALT with effect allele or genomic-reference ALT',
                         'evidence':'author clarification; no new body or independent allele audit','source':AUTHOR},
 'definition_gate_recommendation':'REPORTED_HARE_ONLY_DEFINITION_STOP_RESOLVED_AT_PUBLIC_DOCUMENTED_DEFINITION_LEVEL',
 'hard_gate_caveat':'If frozen gate additionally requires distinct-day semantics, minimum lookback or historical executable phenotyping code, keep those specific requirements unresolved.',
 'admission_scope':'Explicit clinical-insomnia transport; no exact symptom-phenotype substitution or replication claim',
 'family_changes':False,'threshold_changes':False,'new_outcomes_examined':False,'new_gwas_bodies_downloaded':False,
 'native_gwas_compute_performed':False,
 'supplement_pdf':{'url':SUPP,'local_path':str(CACHE/'gia_supplement.pdf'),'bytes':2587868,
                   'md5':digest(CACHE/'gia_supplement.pdf','md5'),'sha256':digest(CACHE/'gia_supplement.pdf'),
                   'official_md5_match':True,'pages_rendered_and_visually_inspected':[2,3]},
 'summary_statistic_license':catalog['terms_of_license'],'article_metadata_license':pmc['license_code'],
 'individual_cohort_overlap':'No participant identifiers audited; no certificate of zero individual overlap inferred'
}
dump('provenance_v2_definition_contract.json', contract)

initial_urls = {
 'basic.pdf':README,'gia_directory.html':'https://ftp.ncbi.nlm.nih.gov/dbgap/studies/phs002453/phs002453.v1.p1/analyses/GIA/',
 'study.html':'https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs002453.v1.p1',
 'manifest.html':DICT,'documents.html':'https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/document.cgi?phd=8760&study_id=phs002453.v1.p1',
 'analyses.html':'https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/analysis.cgi?pha=11340&study_id=phs002453.v1.p1',
 'pmc.html':'https://pmc.ncbi.nlm.nih.gov/articles/PMC12857194/',
 'gia_toc.txt':'https://ftp.ncbi.nlm.nih.gov/dbgap/studies/phs002453/phs002453.v1.p1/analyses/GIA/phs002453.MVP_R4.1000G_AGR.GIA.PheCodes_Neurological.analysis-PI.MULTI.tar.table_of_contents.txt',
 'all_documents.txt':'https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/GetListOfAllObjects.cgi?study_id=phs002453.v1.p1&object_type=document',
 'all_analyses.txt':'https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/GetListOfAllObjects.cgi?study_id=phs002453.v1.p1&object_type=analysis',
 'gia_supplement_response.bin':'https://pmc.ncbi.nlm.nih.gov/articles/instance/12857194/bin/NIHMS2091849-supplement-Supplemental_Text.pdf',
 'publisher_supplement_response.bin':'https://www.science.org/doi/suppl/10.1126/science.adj1182/suppl_file/science.adj1182_sm.pdf',
 'phecodes.html':'https://phewascatalog.org/phecodes', 'exact_gia_analysis.html':GIA,
 'pmc_oa.xml':'https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi?id=PMC12857194',
 'saige_repository.json':'https://api.github.com/repos/exascale-genomics/SAIGE-GPU/git/trees/master?recursive=1',
 'catalog.json':CATALOG,'catalog_metadata.yaml':YAML,
 'uc_repository.html':'https://escholarship.org/uc/item/53c9n629',
 'europe_pmc_search.json':'https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:39024449&format=json',
 'europe_pmc_fulltext.xml':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12857194/fullTextXML',
 'saige_repo_metadata.json':'https://api.github.com/repos/exascale-genomics/SAIGE-GPU',
 'phewas_repo_metadata.json':'https://api.github.com/repos/PheWAS/PheWAS',
 'phecodes_map_response.bin':'https://phewascatalog.org/files/Phecode_map_v1_2_icd10cm_beta.csv.zip',
 'author_response.pdf':AUTHOR,
 'saige_tree.json':'https://api.github.com/repos/exascale-genomics/SAIGE-GPU/git/trees/main?recursive=1',
 'phewas_tree.json':'https://api.github.com/repos/PheWAS/PheWAS/git/trees/master?recursive=1',
 'europe_supplement_response.bin':'https://europepmc.org/articles/PMC12857194/bin/NIHMS2091849-supplement-Supplemental_Text.pdf',
 'gia_supplement_cdn.pdf':'https://cdn.ncbi.nlm.nih.gov/pmc/blobs/324c/12857194/afe2d8e36ab1/NIHMS2091849-supplement-Supplemental_Text.pdf',
 'bioc_list.xml':'https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/supplmat.cgi/bioc_xml/PMC12857194/list',
 'supplement_archive_head.txt':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12857194/supplementaryFiles?inlineImages=false',
 'oa_new.xml':'https://pmc.ncbi.nlm.nih.gov/utils/oa/oa.fcgi?id=PMC12857194',
 'bioc_list_retry.xml':'https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/supplmat.cgi/bioc_xml/39024449/list',
 'europe_supplement_archive_range.bin':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12857194/supplementaryFiles?inlineImages=false',
}
for name in ['phecode_map.rda','phecode_exclude.rda','phecode_rollup_map.rda']:
    initial_urls[name]='https://raw.githubusercontent.com/PheWAS/PheWAS/55dd1c24e228851922400cfba8d7db474565ccc7/data/'+name
initial_urls['phecode_map_docs.txt']='https://raw.githubusercontent.com/PheWAS/PheWAS/55dd1c24e228851922400cfba8d7db474565ccc7/man/phecode_map.Rd'
initial_urls['phecode_create_docs.txt']='https://raw.githubusercontent.com/PheWAS/PheWAS/55dd1c24e228851922400cfba8d7db474565ccc7/man/createPhenotypes.Rd'
ssd_urls={
 'preprint_gia.xml':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10327290/fullTextXML',
 'pmc_cloud_documentation.html':'https://pmc.ncbi.nlm.nih.gov/tools/pmcaws/',
 'pmc_cloud_README.txt':'https://pmc-oa-opendata.s3.amazonaws.com/README.txt',
 'pmc_article_objects.xml':'https://pmc-oa-opendata.s3.amazonaws.com/?list-type=2&prefix=PMC12857194&max-keys=100',
 'gia_supplement.pdf':SUPP,'pmc_article_metadata.json':SUPPMETA,
 'cipher_example.html':'https://phenomics.va.ornl.gov/web/cipher/phenotype-viewer/details?id=15871',
 'cipher_viewer.html':'https://phenomics.va.ornl.gov/web/cipher/phenotype-viewer',
 'cipher_resource_article.xml':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11031216/fullTextXML',
 'cipher_main.js':'https://phenomics.va.ornl.gov/web/main-FFBMV5A5.js',
 'cipher_services.js':'https://phenomics.va.ornl.gov/web/chunk-QPHHAAXS.js',
 'cipher_search.js':'https://phenomics.va.ornl.gov/web/chunk-N4Y3M3KC.js',
 'cipher_phenotype_service.js':'https://phenomics.va.ornl.gov/web/chunk-S66AZSEP.js',
 'cipher_insomnia_search.json':'https://phenomics.va.ornl.gov/web/api/search/phenotypes?query=insomnia&offset=0&limit=20',
 'cipher_example_metadata.json':'https://phenomics.va.ornl.gov/web/api/phenotype/15871',
 'cipher_insomnia_metadata.json':CIPHER,
}
removed={
 'all_analyses.txt':{'bytes':4243123,'sha256':'0b6b94067e75b7186a5641121c036fd8e52333da69fba5e3268f5fe4352ef9de'},
 'manifest.html':{'bytes':1623431,'sha256':'3491a7f2d7c4a4b3d82316a5d84e45b0bc931ea9b7e1dfb3069d996c03ef33d1'}
}
access=[]
for root, urls in [(INITIAL,initial_urls),(CACHE,ssd_urls)]:
    for name,url in urls.items():
        body=root/name; headers=root/(name+'.headers')
        h=headers.read_text(errors='replace') if headers.exists() else ''
        statuses=re.findall(r'^HTTP/\S+\s+(\d+)',h,re.M)
        types=re.findall(r'^content-type:\s*(.*)',h,re.M|re.I)
        dates=re.findall(r'^date:\s*(.*)',h,re.M|re.I)
        record={'resource':name,'url':url,'request_method':'GET','http_status':int(statuses[-1]) if statuses else None,
                'response_content_type':types[-1].strip() if types else None,
                'server_date_header':dates[-1].strip() if dates else None,
                'response_body_retained':body.exists(),'local_cache_path':str(body) if body.exists() else None,
                'response_header_sha256':digest(headers) if headers.exists() else None}
        if body.exists():
            data=body.read_bytes();record.update(bytes=len(data),sha256=hashlib.sha256(data).hexdigest(),
                pdf_magic_valid=data.startswith(b'%PDF-'),snapshot_file_mtime_utc=datetime.fromtimestamp(body.stat().st_mtime,timezone.utc).isoformat())
        elif name in removed:
            record.update(removed[name]);record['cache_disposition']='Redundant own new cache removed during disk exhaustion after SHA256 capture'
        if name=='supplement_archive_head.txt':
            record.update(request_method='HEAD',curl_exit_code=28,request_timeout_seconds=45,bytes=0,
                          access_failure='No response bytes before timeout')
        if name=='europe_supplement_archive_range.bin':
            record.update(request_range='bytes=0-1048575',max_filesize_bytes=1048576,curl_exit_code=28,
                          request_timeout_seconds=30,bytes=0,access_failure='No response bytes before timeout')
        if 'supplement_response' in name or name=='gia_supplement_cdn.pdf':
            record['valid_methods_pdf']=False
        if name=='gia_supplement_response.bin':record['access_failure']='HTTP200 HTML browser challenge; not PDF'
        if name=='gia_supplement_cdn.pdf':record['prior_local_write_failure']='Before retained404 response: curl exit23 / no-space-left-on-device; failed attempt supplied no usable response'
        if name=='cipher_insomnia_search.json':
            record.update(request_method='POST',request_body='[]',request_body_sha256=hashlib.sha256(b'[]').hexdigest(),
                          purpose='Read-only anonymous public phenotype metadata search; no private credentials or cookies')
        if name=='gia_supplement.pdf':record.update(valid_methods_pdf=True,md5=digest(body,'md5'),official_md5_match=True)
        if name=='author_response.pdf':record.update(md5=digest(body,'md5'),valid_pdf=True,pages_visually_inspected=[2])
        access.append(record)

receipt={
 'review_completed_utc':datetime.now(timezone.utc).isoformat(),
 'sealed_v1_commit':'c94ca550875b5ddeaae0b7f6eb99d861c4150d16',
 'scope':'Public metadata and phenotype methods only; no new GWAS body or pair outcomes; no native genetic computation',
 'source_definition_gap_resolution':'Official GIA PDF + official CIPHER exact insomnia mapping/no-exclusion rules + exact GIA EUR count linkage',
 'source_definition_gate':'Resolved at documented public definition level; implementation fields remain explicitly unknown',
 'original_sources_modified':False,'v1_outputs_modified':False,'restricted_access_attempted':False,'contacts_sent':False,
 'verification':{'sealed_provenance_v1_SHA256SUMS':'17 files verified OK on 2026-10-09',
                 'git_diff_sleep_unified_research_v1':'empty on 2026-10-09',
                 'review_v2_checksum_verification':'all generated review files verified OK',
                 'official_icd_mapping_schema':'18 unique code rows; 8 ICD9CM and 10 ICD10CM'},
 'new_summary_statistic_bodies_downloaded':False,'new_result_spreadsheets_downloaded':False,
 'anonymous_public_api_used':True,'challenge_solved':False,
 'artifact_generation':'Metadata serialization plus read-only reference R-data extraction; no estimator, pair tests or thresholds',
 'pdf_inspection':{'gia_methods_pages_rendered':[2,3],'author_response_page_rendered':[2]},
 'disk_incident':'Only redundant own new metadata cache and rendered images were removed; subsequent fetches routed to authorized new SSD namespace',
 'source_cache_root':str(CACHE),'access_records':access,
 'non_http_failures':[{'operation':'createBrowserTab iab for public CIPHER phenotype search','result':'Browser is not available: iab','ui_interaction_performed':False}],
 'official_primary_evidence_files':{
   'gia_methods_pdf':{'path':str(CACHE/'gia_supplement.pdf'),'sha256':digest(CACHE/'gia_supplement.pdf')},
   'cipher_insomnia_record':{'path':str(CACHE/'cipher_insomnia_metadata.json'),'sha256':digest(CACHE/'cipher_insomnia_metadata.json')},
   'gia_analysis':{'path':str(INITIAL/'exact_gia_analysis.html'),'sha256':digest(INITIAL/'exact_gia_analysis.html')},
   'catalog_metadata':{'path':str(INITIAL/'catalog.json'),'sha256':digest(INITIAL/'catalog.json')},
   'author_response':{'path':str(INITIAL/'author_response.pdf'),'sha256':digest(INITIAL/'author_response.pdf')}
 }
}
dump('provenance_v2_source_access_receipt.json',receipt)

paths=sorted(p for p in REVIEW.iterdir() if p.is_file() and p.name!='provenance_v2_SHA256SUMS')
(REVIEW/'provenance_v2_SHA256SUMS').write_text(''.join(f'{digest(p)}  {p.name}\n' for p in paths))
print(json.dumps({'claim_rows':len(rows),'official_mapping_rows':len(mapping),'access_records':len(access),
                  'supplement_md5_verified':True,'files_hashed':len(paths),'new_gwas_outcomes_read':False}))
