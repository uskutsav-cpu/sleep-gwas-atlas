import csv,json,hashlib
from pathlib import Path
base=Path(__file__).resolve().parents[2];here=Path(__file__).parent
paths={'gtex':base/'variant_context/gtex_v8/significant_associations.tsv','predict':base.parent/'brain6-alphagenome-pilot/results/jbrowse_expanded/expanded_variant_evidence.tsv','followup':base/'replication/clinical_three_trait/completed_results.tsv'}
def read(f):
 with f.open(newline='') as h:return list(csv.DictReader(h,delimiter='\t'))
a=read(paths['gtex']); b=read(paths['predict']); c=read(paths['followup']); q=[x for x in a if x['snp']=='rs1476535' and x['type']=='Eqtl']; pred=[x for x in b if x['snp']=='rs1476535']; follow=[x for x in c if x['variant']=='rs1476535'];assert len(q)==3 and len(pred)==len(follow)==1
v=pred[0];g=follow[0];assert v['atlas_variant']==g['grch38']=='chr7:114430980:C>T';assert all(x['variant_id']=='chr7_114430980_C_T_b38' and x['gene_symbol']=='FOXP2' and x['gencode_id'].split('.')[0]=='ENSG00000128573' for x in q)
assert {x['tissue'] for x in q}=={'Pituitary','Adrenal_Gland','Adipose_Subcutaneous'}
rows=[{'variant':v['atlas_variant'],'tissue':x['tissue'],'gene':x['gencode_id'],'gtex_alt_nes':x['normalized_effect_size'],'gtex_nominal_p':x['p_value'],'gtex_response_sha256':x['query_sha256'],'mvp_insomnia_beta_alt':g['insomnia_beta_alt'],'finngen_adhd_beta_alt':g['adhd_beta_alt'],'finngen_hfrs_beta_alt':g['hfrs_beta_alt'],'adult_frontal_dnase_501_log2_alt_ref':v['dnase_log2_ratio_501bp']} for x in q]
with (here/'tissue_crosswalk.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=rows[0],delimiter='\t');w.writeheader();w.writerows(rows)
(here/'receipt.json').write_text(json.dumps({'inputs_sha256':{k:hashlib.sha256(v.read_bytes()).hexdigest() for k,v in paths.items()},'matched_tissues':len(rows)},indent=2)+'\n')
for r in rows:print(r)
