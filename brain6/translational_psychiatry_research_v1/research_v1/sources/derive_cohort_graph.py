#!/usr/bin/env python3
"""Named cohort graph, explicit unknown pair overlap, two-source ADHD checks."""
import csv,hashlib,json
from pathlib import Path
import openpyxl
ROOT=Path(__file__).resolve().parents[4];AREA=ROOT/'brain6/translational_psychiatry_research_v1/research_v1';RAW=ROOT/'work/tp-literature-20261007'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,rr):
 with p.open('w',newline='') as f:
  w=csv.DictWriter(f,list(rr[0]),delimiter='\t');w.writeheader();w.writerows(rr)
def main():
 a=RAW/'adhd_cohort_summary.xlsx';b=RAW/'demontis_2023_supp_tables.xlsx'
 w=openpyxl.load_workbook(a,read_only=True,data_only=True)
 fa={str(r[0]).strip():(i,r) for i,r in enumerate(w['cohortData'].values,1) if isinstance(r[1],int) and isinstance(r[2],int)}
 w=openpyxl.load_workbook(b,read_only=True,data_only=True);rr=[];edges=[];checks=[]
 for i,r in enumerate(w['1.Cohorts included'].values,1):
  if not isinstance(r[1],int) or not isinstance(r[3],int):continue
  name=r[0].strip();fi,fr=fa[name];assert fr[1]==r[1] and fr[2]==r[3]
  rr.append(dict(cohort=name,cases=r[1],controls=r[3],female_fraction_cases=r[2],female_fraction_controls=r[4],age_group=r[5],sample_design=r[6],ancestry=r[7],covariates=r[8],source_table='Demontis2023 Supplementary Table1',source_row=i,source_sha256=sha(b),figshare_sheet='cohortData',figshare_row=fi,figshare_sha256=sha(a),cross_source_count_match=True,overlap_with_FinnGen='UNVERIFIED; named roster is not a participant-level or control-pool linkage'))
  edges.append(dict(from_node='Brain6_ADHD_Demontis2023',to_node='cohort:'+name,relation='CONTRIBUTING_COHORT',evidence_level='SOURCE_ATTESTED',cases=r[1],controls=r[3],source=str(b.relative_to(ROOT)),source_sha256=sha(b),limitation='Clinical case/control membership and secondary cohort aliases require source-level participant comparison'))
 assert len(rr)==13 and sum(x['cases'] for x in rr)==38691 and sum(x['controls'] for x in rr)==186843
 checks.append(dict(check='Two primary-source cohort manifests agree on all13 named cohorts and case/control counts',observed='13/13 names and26/26 counts match; summed cases38691 controls186843',status='PASS',source_sha256=sha(a)+';'+sha(b)))
 for node,cohort in [('Brain6_INS_Jansen2019_public','UK_Biobank'),('Brain6_long_sleep_Dashti2019','UK_Biobank'),('Brain6_MDD_Howard2019_public','UK_Biobank'),('Brain6_MDD_Howard2019_public','PGC_139k')]:
  p=ROOT/'brain6/confirmatory_v2/SOURCE_AND_ACCESS_LEDGER.tsv'
  edges.append(dict(from_node=node,to_node='cohort:'+cohort,relation='CONTRIBUTING_COHORT',evidence_level='HISTORICAL_SOURCE_LEDGER_ATTESTED',cases='',controls='',source=str(p.relative_to(ROOT)),source_sha256=sha(p),limitation='Named component establishes documented UKB overlap; individual counts/covariance still unknown'))
 biobanks=['Arctic Biobank','Auria Biobank','Biobank Borealis of Northern Finland','Biobank of Eastern Finland','Central Finland Biobank','Finnish Red Cross Blood Service Biobank','Finnish Clinical Biobank Tampere','Helsinki Biobank','Terveystalo Biobank','THL Biobank']
 p=RAW/'finngen_cohorts.html'
 for name in biobanks:
  assert name in p.read_text()
  edges.append(dict(from_node='FinnGen_R13',to_node='organisation:'+name,relation='PARTICIPATING_BIOBANK',evidence_level='SOURCE_ATTESTED_ORGANISATION_NOT_ENDPOINT_SAMPLE_MANIFEST',cases='',controls='',source=str(p.relative_to(ROOT)),source_sha256=sha(p),limitation='Organisation and cohort names are different units; endpoint-specific participants and legacy controls unlinked'))
 for code in ['F5_INSOMNIA','F5_ADHD','F5_DEPRESSIO']:
  p=RAW/'finngen_r13_phenotypes.json'
  edges.append(dict(from_node='FinnGen_R13_'+code,to_node='FinnGen_R13',relation='ENDPOINT_SUBSET',evidence_level='SOURCE_ATTESTED',cases='',controls='',source=str(p.relative_to(ROOT)),source_sha256=sha(p),limitation='Cases/comorbidities and controls overlap across endpoints; exact pairwise overlap not provided'))
 vfile=ROOT/'work/tp-finngen-fixed-20261007/rs77960.json';annot=json.loads(vfile.read_text())['variant']['annotation']['annot']
 for token in ['NFBC66','NFBC86']:
  assert any(token in k for k in annot)
  edges.append(dict(from_node='FinnGen_public_variant_annotation',to_node='legacy_label:'+token,relation='LEGACY_GENOTYPE_ANNOTATION_LABEL',evidence_level='SOURCE_ATTESTED_LABEL_ONLY',cases='',controls='',source=str(vfile.relative_to(ROOT)),source_sha256=sha(vfile),limitation='Annotation AF/INFO dataset label is not proof of endpoint case/control membership; no matching ADHD cohort name in current roster does not establish participant disjointness'))
 unknown=[]
 for d in ['Brain6_ADHD_Demontis2023','Brain6_MDD_Howard2019_public','Brain6_INS_Jansen2019_public']:
  unknown.append(dict(discovery_source=d,external_source='FinnGen_R13',overlap_classification='UNVERIFIED',named_roster_result='No exact ADHD cohort name matches ten declared FinnGen biobank organisations; different naming units preclude disjointness conclusion' if 'ADHD' in d else 'No complete endpoint-specific participant manifest links discovery and FinnGen',case_control_overlap='UNVERIFIED',reason='No person-level linkage, control-pool manifest, complete legacy aliases or endpoint-specific cohort allocation available',allowed_class='PHENOTYPE_SENSITIVITY',next_required_evidence='Consortium confirmation or de-identified source-level overlap counts and complete recruitment/control-pool alias mapping'))
 save(AREA/'sources/ADHD_DISCOVERY_COHORTS.tsv',rr);save(AREA/'sources/SOURCE_COHORT_GRAPH.tsv',edges);save(AREA/'sources/COHORT_OVERLAP_DECISIONS.tsv',unknown);save(AREA/'sources/cohort_cross_source_checks.tsv',checks)
 print('ADHD13 cohorts, tenPGC names, sums38691/186843 verified; graph edges',len(edges),'overlap remains UNVERIFIED')
if __name__=='__main__':main()
