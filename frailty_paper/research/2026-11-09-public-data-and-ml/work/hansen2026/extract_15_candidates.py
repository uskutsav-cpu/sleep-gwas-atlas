import csv,gzip,hashlib,json
from pathlib import Path
p=Path('work/hansen2026/parallel_verified.tsv.gz');h=hashlib.md5()
with p.open('rb') as f:
 for b in iter(lambda:f.read(1048576),b''):h.update(b)
assert h.hexdigest()=='f634c1ceadf8bd32e98b311aab8a973f'
rows=list(csv.DictReader(open('outputs/integrated_15_candidate_review.csv')));keys={(r['GRCh37_chr'],r['GRCh38_pos']):r for r in rows};found={k:[] for k in keys};n=0
with gzip.open(p,'rt') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  n+=1;k=(r['chromosome'],r['base_pair_location'])
  if k in found:found[k].append(r)
out=[]
for k,c in keys.items():
 a=c['counted_allele'].upper();b=c['other_allele'].upper();match=[r for r in found[k] if {r['effect_allele'].upper(),r['other_allele'].upper()}=={a,b}];row={'variant':c['variant'],'chromosome':k[0],'position_grch38':int(k[1]),'counted_allele':a,'other_allele':b,'position_rows':len(found[k]),'exact_matches':len(match)}
 if len(match)==1:
  x=match[0];sign=1 if x['effect_allele'].upper()==a else -1;ref=x['effect_allele'] if x['ref_allele']=='EA' else x['other_allele'] if x['ref_allele']=='OA' else None;bits=x['variant_id'].split('_');row.update({'beta_counted_unverified_units':sign*float(x['beta']),'se_unverified_units':float(x['standard_error']),'p':float(x['p_value']),'source_effect_allele':x['effect_allele'],'source_eaf':x['effect_allele_frequency'],'source_variant_id':x['variant_id'],'source_ref':ref,'source_ref_agrees_variant_id':len(bits)==4 and ref==bits[2],'sample_size':x['n']})
 out.append(row)
Path('outputs/stage6_hansen_15_candidate_screen.json').write_text(json.dumps({'source_md5':h.hexdigest(),'source_records':n,'rows':out},indent=2)+'\n');print('rows',n,'matched',sum(r['exact_matches']==1 for r in out),'p<.05',[(r['variant'],r['p']) for r in out if 'p' in r and r['p']<.05])
