#!/usr/bin/env python3
"""Exact primary supplement extraction; coordinate overlap is not causal identity.

Requires openpyxl (bundled Codex Python runtime). No outcomes downloaded here.
Worksheet and row references are one based. Candidate windows come from frozen
pair-specific IDs/crosswalk, not geographic merged display intervals.
"""
import csv,hashlib,json,re,zipfile,xml.etree.ElementTree as E
from pathlib import Path
import openpyxl
ROOT=Path(__file__).resolve().parents[4]
AREA=ROOT/'brain6/translational_psychiatry_research_v1/research_v1';RAW=ROOT/'work/tp-literature-20261007'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,rows,fields=None):
 if fields is None:fields=list(rows[0])
 with p.open('w',newline='') as f:
  w=csv.DictWriter(f,fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(rows)
def rows_at(w,s,h):
 rr=list(w[s].values);headers=rr[h-1]
 return [(i,{str(k):v for k,v in zip(headers,row) if k is not None}) for i,row in enumerate(rr[h:],h+1)]
def main():
 p=RAW/'lin_supp1.xlsx';w=openpyxl.load_workbook(p,read_only=True,data_only=True)
 t8={(d['Trait_pair'],d['rsID']):(i,d) for i,d in rows_at(w,'Supplementary Table 8',2)}
 t7={(d.get('Trait pair'),d.get('Top SNP')):(i,d) for i,d in rows_at(w,'Supplementary Table 7',3)}
 signals=[]
 for i,d in rows_at(w,'Supplementary Table 6',2):
  if d['Trait_pair'] not in ['INS-ADHD','INS-MDD']:continue
  rr,cc=t8.get((d['Trait_pair'],d['rsID']),(None,{}));br,bb=t7.get((d['Trait_pair'],d['rsID']),(None,{}))
  signals.append(dict(study='Lin_2026',doi='10.1007/s00335-026-10232-5',pair_id='insomnia__'+d['Trait_pair'].split('-')[1].lower(),rsid=d['rsID'],chrom=int(d['CHR']),pos=int(d['BP']),start=int(d['start']),end=int(d['end']),build='GRCh37',build_evidence='Public Supplementary Methods FUMA v1.5.2 explicitly uses GRCh37',allele1=d['A1'],allele2=d['A2'],placo_p=d['P.PLACO'],nearest_gene=d['Nearest gene'],source_sheet='Supplementary Table 6',source_row=i,coloc_sheet='Supplementary Table 8' if cc else '',coloc_row=rr,h3=cc.get('PP.H3.abf'),h4=cc.get('PP.H4.abf'),coloc_nsnp=cc.get('nsnps'),coloc_top_snp=cc.get('Top.causal.SNP'),beta_sleep=bb.get('Beta_Sleep'),beta_outcome=bb.get('Beta_PSY'),p_sleep=bb.get('P_Sleep'),p_outcome=bb.get('P_PSY'),effect_sheet='Supplementary Table 7',effect_row=br,source_sha256=sha(p),interpretation='Published same trait-pair signal; single-causal ABF is descriptive, not independent validation'))
 ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
 p=RAW/'jia_supp1.docx';root=E.fromstring(zipfile.ZipFile(p).read('word/document.xml'))
 tables=[]
 for table in root.findall('.//w:tbl',ns):
  tables.append([[''.join(t.text or '' for t in c.findall('.//w:t',ns)) for c in row.findall('w:tc',ns)] for row in table.findall('w:tr',ns)])
 trait_map={'Insomnia-MDD':'insomnia__mdd','Long sleep duration-BP':'longsleep__bipolar','Long sleep duration-SCZ':'longsleep__scz'}
 coloc={(row[0],row[1]):(i,row) for i,row in enumerate(tables[6],1) if len(row)>=6 and row[0] in trait_map}
 for i,d in enumerate(tables[2],1):
  if len(d)<11 or d[1] not in trait_map:continue
  ci,c=coloc.get((d[1],d[2]),(None,[]))
  signals.append(dict(study='Jia_2025',doi='10.1093/sleep/zsae209',pair_id=trait_map[d[1]],rsid=d[2],chrom=int(d[3]),pos=int(d[4]),start=None,end=None,build='GRCh37_COORDINATE_CONSISTENT',build_evidence='Source does not expressly attest build; selected lead Ensembl37 positions verified separately',allele1='',allele2='',placo_p=float(d[7]),nearest_gene=d[5],source_sheet='Supplementary Table S3',source_row=i,coloc_sheet='Supplementary Table S7',coloc_row=ci,h3=None,h4=float(c[3]) if c else None,coloc_nsnp=None,coloc_top_snp=c[4] if c else None,beta_sleep=None,beta_outcome=None,p_sleep=None,p_outcome=None,effect_sheet='',effect_row=None,source_sha256=sha(p),interpretation='Alleles unreported; no exact biallelic identity or shared-causal inference from rsID/proximity alone; long sleep is binary Dashti PMID30846698 per S1, exact Brain6 source harmonization separate'))
 # Zu C is an exact HTML table fact already independently acquired in the old cache.
 zp=ROOT/'work/literature-20261007/zu_table3.html'
 signals.append(dict(study='Zu_2026',doi='10.1038/s41398-026-04166-4',pair_id='insomnia__adhd',rsid='rs77960',chrom=5,pos=103964585,start=None,end=None,build='GRCh37_COORDINATE_VERIFIED',build_evidence='Official Ensembl37 mapping matches published Table 3 coordinate; full source build provenance must not be inferred from one rsID',allele1='',allele2='',placo_p=6.54e-15,nearest_gene='RP11-6N13.1',source_sheet='Main Table 3',source_row=14,coloc_sheet='Main Table 3',coloc_row=14,h3=.01,h4=.99,coloc_nsnp=None,coloc_top_snp=None,beta_sleep=None,beta_outcome=None,p_sleep=3.34e-8,p_outcome=2.46e-13,effect_sheet='',effect_row=None,source_sha256=sha(zp),interpretation='Reported beta_ADHD .9294 resembles OR; orientation not certified; current EUR reference LD comparison is separate evidence'))
 save(AREA/'novelty/published_candidate_leads.tsv',signals)
 # Schipper has side-by-side tables. Explicit offsets prevent wrong-window parsing.
 p=RAW/'schipper_v2_supp_tables.xlsx';w=openpyxl.load_workbook(p,read_only=True,data_only=True)
 sch=[]
 for i,d in enumerate(w['6. LAVA'].values,1):
  if d[1:4]!=(11,112755447,113889019):continue
  sch.append(dict(study='Schipper_v2_20260311',doi='10.1101/2025.10.18.25338281',source_sheet='6. LAVA',source_row=i,pair='INS_DEP',chrom=11,start=112755447,end=113889019,n_snps=d[4],n_pcs=d[5],local_rg=d[6],rg_lower=d[7],rg_upper=d[8],p=d[12],h0=None,h1=None,h2=None,h3=None,h4=None,n_coloc_snps=None,source_sha256=sha(p)))
 for i,d in enumerate(w['7. Coloc'].values,1):
  if d[24:27]!=(11,112755447,113889019):continue
  sch.append(dict(study='Schipper_v2_20260311',doi='10.1101/2025.10.18.25338281',source_sheet='7. Coloc right table',source_row=i,pair='INS_DEP' if d[16]=='LAVA_EUR_INS' else 'DEP_ANX',chrom=11,start=112755447,end=113889019,n_snps=None,n_pcs=None,local_rg=None,rg_lower=None,rg_upper=None,p=None,h0=d[17],h1=d[18],h2=d[19],h3=d[20],h4=d[21],n_coloc_snps=d[22],source_sha256=sha(p)))
 save(AREA/'novelty/schipper_locus_evidence.tsv',sch)
 old=list(csv.DictReader((ROOT/'brain6/confirmatory_v2/review/locus_novelty_crosswalk.tsv').open(),delimiter='\t'))
 planned=[dict(candidate_locus_id='candidate_A_INS_ADHD',pair_id='insomnia__adhd',interval_grch37='chr11:112460000-114260000',Brain6_lead=''),dict(candidate_locus_id='candidate_A_INS_MDD',pair_id='insomnia__mdd',interval_grch37='chr11:112460000-114260000',Brain6_lead=''),dict(candidate_locus_id='candidate_B_INS_MDD',pair_id='insomnia__mdd',interval_grch37='chr6:100630000-102640000',Brain6_lead='')]
 master=[]
 for c in old+planned:
  chrom,start,end=map(int,re.match(r'chr(\d+):(\d+)-(\d+)',c['interval_grch37']).groups())
  match=[s for s in signals if s['pair_id']==c['pair_id'] and s['chrom']==chrom and start<=s['pos']<=end]
  has_sch= c['pair_id']=='insomnia__mdd' and chrom==11 and start<=113889019 and end>=112755447
  common=dict(candidate_locus_id=c['candidate_locus_id'],pair_id=c['pair_id'],pair_window_grch37=c['interval_grch37'],brain6_lead=c['Brain6_lead'],competing_studies=';'.join(sorted({s['study'] for s in match}|({'Schipper_v2_20260311'} if has_sch else set()))),published_leads=';'.join(s['rsid'] for s in match),exact_rsid_recurrence=any(s['rsid'] in c['Brain6_lead'].split(';') for s in match),source_tables=';'.join(f"{s['study']}:{s['source_sheet']}:row{s['source_row']}" for s in match),source_sha256=';'.join(sorted({s['source_sha256'] for s in match}|({sch[0]['source_sha256']} if has_sch else set()))),classification='KNOWN_SAME_PAIR_REGION' if match or has_sch else 'UNRESOLVED_NO_NOVELTY_INFERENCE',distinction='Prior signal at same pair/window; causal identity requires allele/LD/multisignal analysis' if match else ('Prior INS_DEP local correlation at overlapping window; not equivalent to shared-causal support' if has_sch else 'This targeted review does not establish novelty; long-sleep phenotype and other loci need dedicated supplement comparisons'),allowed_wording='Previously reported region' if match or has_sch else 'Candidate region with novelty unresolved',prohibited_wording='Novel locus; independently replicated causal variant; nearest gene is causal',review_cutoff='2026-10-07')
  master.append(common)
 save(AREA/'NOVELTY_MASTER.tsv',master)
 # Build/allele verification does not require whole-genome downloads.
 mappings=[]
 for s in signals:
  p=RAW/(s['rsid']+'_grch37.json')
  if not p.exists():p=ROOT/'work/literature-20261007'/p.name
  if not p.exists():continue
  d=json.loads(p.read_text());mm=[m for m in d['mappings'] if m['assembly_name']=='GRCh37' and m['seq_region_name']==str(s['chrom'])]
  if not mm:continue
  m=mm[0];mappings.append(dict(study=s['study'],rsid=s['rsid'],published_pos=s['pos'],ensembl_pos=m['start'],coordinate_match=m['start']==s['pos'],strand=m['strand'],ensembl_alleles=m['allele_string'],published_allele1=s['allele1'],published_allele2=s['allele2'],source_sha256=sha(p),interpretation='Coordinate consistency only; complement/order and source effect definition separate'))
 save(AREA/'novelty/source_build_checks.tsv',mappings)
 (AREA/'novelty/extraction_environment.json').write_text(json.dumps(dict(openpyxl_version=openpyxl.__version__,script_sha256=sha(Path(__file__)),signal_rows=len(signals),protected_candidate_rows=len(old),master_rows=len(master)),indent=2)+'\n')
 print('signals',len(signals),'master',len(master),'known',sum(x['classification']=='KNOWN_SAME_PAIR_REGION' for x in master),'build_checks',len(mappings))
if __name__=='__main__':main()
