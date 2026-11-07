#!/usr/bin/env python3
import csv,hashlib,json,math
from pathlib import Path
p=Path(__file__).resolve().parent;atlas=p.parents[1];inputs=Path('/Volumes/Extreme SSD/brain6-work/exploratory-finemap-coloc-v2/inputs')
def read(path):
 with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
prov=json.loads((atlas.parent/'brain6-alphagenome-pilot/results/provenance.json').read_text())['input_sha256'];matches=read(p/'matched_candidates.tsv')
assert {x['rsid'] for x in matches}=={'rs1093463','rs1989903'}
regions={x['rsid']:x for x in read(atlas/'variant_context/crossdisorder_mpra/ranked_coverage.tsv')};out=[];source={}
for x in matches:
 rs=x['rsid'];meta=regions[rs];chrom,posref=meta['source_key'].split(':',1);pos,alleles=posref.split('_');ref,alt=alleles.split(':')
 assert chrom=='chr'+x['chromosome'] and {ref,alt}=={x['allele1'],x['allele0']}
 file=inputs/(x['ranked_candidate_region']+'__insomnia.tsv');actual=sha(file);assert actual==prov[file.name];source[file.name]=actual
 matches_source=[r for r in read(file) if r['SNP']==rs];assert len(matches_source)==1
 row=matches_source[0];assert row['CHR']==x['chromosome'] and row['POS']==pos and {row['A1'],row['A2']}=={ref,alt}
 db=float(row['BETA'])*(1 if row['A1']==alt else -1);ds=float(row['SE']);eb=float(x['beta_allele1'])*(1 if x['allele1']==alt else -1);es=float(x['se'])
 assert all(math.isfinite(v) for v in (db,ds,eb,es)) and ds>0 and es>0
 out.append({'variant':rs,'alt':alt,'discovery_beta_alt':db,'discovery_se':ds,'discovery_ci_low':db-1.96*ds,'discovery_ci_high':db+1.96*ds,'estonian_beta_alt':eb,'estonian_se':es,'estonian_ci_low':eb-1.96*es,'estonian_ci_high':eb+1.96*es,'estonian_p':float(x['p']),'same_sign':db*eb>0})
with (p/'direction_comparison.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(out[0]),delimiter='\t');w.writeheader();w.writerows(out)
r={'source_hashes':source,'input_estonian_sha256':sha(p/'matched_candidates.tsv'),'output_sha256':sha(p/'direction_comparison.tsv'),'rows':out,'limitations':'Both variants and sources post-selected; weak Estonian effects do not establish sign, replication or shared mechanism. Wald intervals ignore selection and phenotype differences; independence unverified.'}
(p/'direction_receipt.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
