import csv,gzip,hashlib,json,math,os
from pathlib import Path
p=Path(os.environ.get('HANSEN_GWAS_PATH','work/hansen2026/GCST90854255.tsv.gz'));  h=hashlib.md5()
with p.open('rb') as f:
 for b in iter(lambda:f.read(1048576),b''):h.update(b)
assert h.hexdigest()=='f634c1ceadf8bd32e98b311aab8a973f','download incomplete or checksum mismatch'
c=list(csv.DictReader(open('outputs/expanded_22_candidate_risk_gene_table.csv')));pos={(r['GRCh38_chr'],r['GRCh38_pos']):r for r in c};hits={k:[] for k in pos};n=0
with gzip.open(p,'rt') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  n+=1;k=(r['chromosome'],r['base_pair_location'])
  if k in hits:hits[k].append(r)
freq={r['variant']:r for r in csv.DictReader(open('outputs/22_score_reference_group_allele_frequencies.csv'))}
reference={r['variant']:r['reference_base'] for r in csv.DictReader(open('outputs/stage6_hansen_candidate_grch38_reference.csv'))}
out=[]
for k,record in pos.items():
 a=record['counted_allele'].upper();b=record['other_allele'].upper();same=[r for r in hits[k] if {r['effect_allele'].upper(),r['other_allele'].upper()}=={a,b}];row={'variant':record['variant'],'chromosome':k[0],'position38':k[1],'candidate_allele':a,'other':b,'source_rows_at_position':len(hits[k]),'allele_pair_matches':len(same),'status':'missing' if not hits[k] else 'allele_mismatch' if not same else 'duplicate' if len(same)>1 else 'matched'}
 if len(same)==1:
  x=same[0]; sign=1 if x['effect_allele'].upper()==a else -1;beta=sign*float(x['beta']); source_eaf=float(x['effect_allele_frequency']); candidate_af=source_eaf if sign==1 else 1-source_eaf; reference_af=float(freq[record['variant']]['reference_AFR_frequency']); allele_fields=x['variant_id'].split('_'); source_ref=x['effect_allele'] if x['ref_allele']=='EA' else x['other_allele'] if x['ref_allele']=='OA' else ''; source_ref_agrees=(len(allele_fields)==4 and source_ref==allele_fields[2]); row.update({'source_effect_allele':x['effect_allele'],'source_other_allele':x['other_allele'],'source_eaf':x['effect_allele_frequency'],'source_ref_allele':x['ref_allele'],'source_variant_id':x['variant_id'],'source_ref_matches_variant_id':source_ref_agrees,'ucsc_grch38_reference_base':reference[record['variant']],'source_ref_matches_ucsc':source_ref==reference[record['variant']],'candidate_af':candidate_af,'reference_AFR_candidate_af':reference_af,'absolute_af_difference':abs(candidate_af-reference_af),'aligned_beta_source_units_unverified':beta,'se_source_units_unverified':x['standard_error'],'p_value':x['p_value'],'n':x['n'],'palindromic':{a,b} in [{'A','T'},{'C','G'}]})
 out.append(row)
fields=list(dict.fromkeys(k for r in out for k in r));
with open('outputs/stage6_hansen_grip_candidate_lookup.csv','w') as f:
 w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(out)
Path('outputs/stage6_hansen_grip_receipt.json').write_text(json.dumps({'source_md5':h.hexdigest(),'source_records':n,'candidate_rows':len(out),'matched':sum(r['status']=='matched' for r in out),'source_url':'https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90854001-GCST90855000/GCST90854255/GCST90854255.tsv.gz','source_not_harmonized':True},indent=2)+'\n')
print(n,{s:sum(r['status']==s for r in out) for s in ['matched','missing','allele_mismatch','duplicate']})
