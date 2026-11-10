import csv,hashlib,json,openpyxl
from pathlib import Path
p=Path('outputs');f=Path('work/clhls_genetic_supp/ACEL-20-e13323-s008.xlsx');s=openpyxl.load_workbook(f,read_only=True,data_only=True).active;it=s.iter_rows(min_row=3,values_only=True);hdr=next(it);assert hdr[1]=='snp' and hdr[6]=='alleles_fwd';source={};n=0
for row in it:
 if not isinstance(row[1],str) or not row[1].startswith('rs'):continue
 n+=1;source.setdefault(row[1],[]).append(row)
lists={}
for tag,fn in [('score15','unvalidated_15_variant_score_prototype.csv'),('score22','unvalidated_22_shared_direction_score.csv')]:
 for r in csv.DictReader((p/fn).open()):lists.setdefault(r['variant'],[]).append((tag,r))
rows=[]
for rs,variants in lists.items():
 found=source.get(rs,[])
 for tag,r in variants:
  x=found[0] if len(found)==1 else None; forward=str(x[6]) if x else '';pair={r['counted_allele'],r['other_allele']};gen=set(forward.split('/'));rows.append({'list':tag,'variant':rs,'list_alleles':'/'.join([r['counted_allele'],r['other_allele']]),'chip_rsid_rows':len(found),'chip_alleles_fwd':forward,'chip_chr':x[4] if x else '', 'chip_position_source_build_unverified':x[5] if x else '', 'exact_unordered_pair':bool(x and len(gen)==2 and gen==pair),'interpretation':'Design inventory only; chip build, strand, sample genotype QC and phenotype linkage remain unverified.'})
with (p/'stage6_clhls_chip_candidate_coverage.csv').open('w',newline='') as h:
 w=csv.DictWriter(h,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
report={'source_url':'https://pmc.ncbi.nlm.nih.gov/articles/PMC7963337/','table':'Supplementary Table S3; designed chip inventory, not QC-passing genotypes','sha256':hashlib.sha256(f.read_bytes()).hexdigest(),'chip_rsid_rows':n,'results':{tag:{'candidates':sum(r['list']==tag for r in rows),'unique_rsid_matches':sum(r['list']==tag and r['chip_rsid_rows']==1 for r in rows),'exact_forward_allele_pair_matches':sum(r['list']==tag and r['exact_unordered_pair'] for r in rows)} for tag in ['score15','score22']},'warning':'Not individual-level score validation; observed or imputed genotypes and build/strand are not verified.'}
(p/'stage6_clhls_chip_candidate_coverage.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
