#!/usr/bin/env python3
import collections,csv,hashlib,json,re,zipfile
from pathlib import Path
import openpyxl
p=Path(__file__).resolve().parent;atlas=p.parents[1];xlsx=p/'media-2.xlsx';archive=p/'supplementary.zip'
with zipfile.ZipFile(archive) as z:assert z.testzip() is None and z.read('media-2.xlsx')==xlsx.read_bytes()
assert hashlib.md5(xlsx.read_bytes()).hexdigest()=='71e7719dd7e4c25de11c151111e1582f'
with zipfile.ZipFile(xlsx) as z:assert z.testzip() is None
with (atlas/'prediction_rank_comparison/variant_ranks.tsv').open() as f:candidates=list(csv.DictReader(f,delimiter='\t'))
with (p.parent/'six_lead_joint/six_leads.tsv').open() as f:leads=list(csv.DictReader(f,delimiter='\t'))
assert len(candidates)==425 and len(leads)==6
rank=collections.defaultdict(list)
for x in candidates:
 m=re.fullmatch(r'chr(\d+):(\d+):([ACGT]+)>([ACGT]+)',x['grch38']);assert m
 ch,pos,ref,alt=m.groups();rank[(int(ch),int(pos),frozenset((ref,alt)))].append(x)
w=openpyxl.load_workbook(xlsx,read_only=True,data_only=True);s=w['ST5'];it=s.iter_rows(values_only=True)
title=next(it);head=next(it);assert 'Estonian Biobank' in title[0] and head[:6]==('Phenotype','CHR','POS','SNPID','ALLELE1','ALLELE0')
rows=[x for x in it if x and x[0]=='Insomnia'];counts=collections.Counter();matched=[];lead_keys={(int(x['chr']),int(x['position38']),x['variant']) for x in leads}
for r in rows:
 ch=int(r[1]);pos=int(r[2]);rs=str(r[3]);a1=str(r[4]);a0=str(r[5]);key=(ch,pos,frozenset((a1,a0)));values=rank.get(key,[])
 if len(values)>1:counts['duplicate_candidate_key']+=1
 if not values:continue
 for v in values:
  counts['matched_candidate_rows']+=1
  matched.append({'rsid':rs,'candidate_rsid':v['snp'],'chromosome':ch,'position38':pos,'allele1':a1,'allele0':a0,'beta_allele1':r[10],'se':r[11],'p':r[14],'n':r[8],'ranked_candidate_region':v['region'],'id_agrees':rs==v['snp'],'is_one_of_six_leads':(ch,pos,rs) in lead_keys})
fields=['rsid','candidate_rsid','chromosome','position38','allele1','allele0','beta_allele1','se','p','n','ranked_candidate_region','id_agrees','is_one_of_six_leads']
with (p/'matched_candidates.tsv').open('w',newline='') as f:writer=csv.DictWriter(f,fields,delimiter='\t');writer.writeheader();writer.writerows(matched)
receipt={'source':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12668095/supplementaryFiles','paper':'https://doi.org/10.1101/2025.11.21.25340447','archive_sha256':hashlib.sha256(archive.read_bytes()).hexdigest(),'workbook_sha256':hashlib.sha256(xlsx.read_bytes()).hexdigest(),'workbook_md5':hashlib.md5(xlsx.read_bytes()).hexdigest(),'europe_pmc_xml_md5_match':True,'archive_and_workbook_zip_integrity':'passed','table':'ST5 Estonian Biobank clinical sleep disorder selected follow-up','insomnia_data_rows':len(rows),'candidate_count':len(candidates),'lead_count':len(leads),'matched':matched,'counts':dict(counts),'rows_sha256':hashlib.sha256((p/'matched_candidates.tsv').read_bytes()).hexdigest(),'limitations':'Preprint, selected lead-only table and post-hoc candidate intersection, not full regional or genome-wide results. No independent ADHD, locus colocalization, or ancestry/sample overlap validation. Effect direction comparison not conducted.'}
(p/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
