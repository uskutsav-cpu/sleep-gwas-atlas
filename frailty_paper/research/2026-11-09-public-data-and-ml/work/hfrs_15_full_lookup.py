import csv,hashlib,collections,json
from pathlib import Path
root=Path('/Volumes/Extreme SSD/sleep-frailty-full-autonomous-v5-20261009/frailty_paper/full_autonomous_research_v5/genetics/new_public_validation/raw')
p=root/'GCST90707291.tsv';expected='38d9434f50a4493def9d2c9d1cc1b4b3'
c={r['variant']:r for r in csv.DictReader(open('outputs/candidate_risk_gene_review_15.csv'))}
b={r['variant']:r for r in csv.DictReader(open('outputs/prototype_variants_build_reference.csv'))}
print('checking',len(c),'candidate markers',flush=True)
found=collections.defaultdict(list);h=hashlib.md5();n=0
with p.open('rb') as f:
 header=None
 for line in f:
  h.update(line)
  if header is None:header=line.decode().rstrip('\n\r').split('\t');continue
  n+=1
  vals=line.decode().rstrip('\n\r').split('\t');rs=vals[-1]
  if rs in c:found[rs].append(dict(zip(header,vals)))
 if h.hexdigest()!=expected:raise ValueError('MD5 mismatch '+h.hexdigest())
out=[]
for rs,r in c.items():
 target=b[rs]; match=[x for x in found[rs] if x['chromosome']==target['grch38_chr'] and x['base_pair_location']==target['grch38_pos']]
 # Determine eligible exact allele contrast; never auto-complement or flip multiallelic variants.
 for x in match:
  a=r['counted_allele'];o=r['other_allele'];ea=x['effect_allele'];oa=x['other_allele']
  if (ea,oa)==(a,o):align='direct'; beta=float(x['beta'])
  elif (ea,oa)==(o,a):align='swap'; beta=-float(x['beta'])
  else:align='mismatch';beta=None
  out.append({'variant':rs,'GRCh38_chr':x['chromosome'],'GRCh38_pos':x['base_pair_location'],'candidate_counted_allele':a,'candidate_other_allele':o,'hfrs_source_effect_allele':ea,'hfrs_source_other_allele':oa,'alignment':align,'hfrs_beta_on_counted_allele':beta,'hfrs_se':x['standard_error'],'hfrs_p':x['p_value'],'hfrs_source_effect_allele_frequency':x['effect_allele_frequency'],'candidate_FI_beta':r['FI_meta_beta'],'candidate_FFS_beta_inferred':r['FFS_beta_inferred_A1'],'candidate_insomnia_UKB_beta':r['insomnia_Jansen_beta'],'rsid_source_rows':len(found[rs]),'position_source_rows':len(match),'reference_multiallelic':r['multiallelic_reference'],'status':'Exploratory selected variants; source allele/build string match only; no score validation or causal inference'})
 if not match:out.append({'variant':rs,'GRCh38_chr':target['grch38_chr'],'GRCh38_pos':target['grch38_pos'],'candidate_counted_allele':r['counted_allele'],'candidate_other_allele':r['other_allele'],'alignment':'no_coordinate_match','rsid_source_rows':len(found[rs]),'position_source_rows':0,'status':'Not found at mapped position; not a null association'})
with open('outputs/candidate_finngen_hfrs_full_gwas_lookup.csv','w',newline='') as f:
 keys=list(dict.fromkeys(k for row in out for k in row));w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(out)
print('rows_scanned',n,'md5',h.hexdigest(),'candidates_with_id',len(found),'aligned',sum(r['alignment'] in ['direct','swap'] for r in out),flush=True)
