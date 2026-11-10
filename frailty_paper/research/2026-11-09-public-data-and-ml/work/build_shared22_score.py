import csv,pathlib,hashlib,json
source=list(csv.DictReader(open('outputs/22_candidate_risk_direction_table.csv')));orig={r['variant']:r for r in csv.DictReader(open('outputs/unvalidated_22_variant_score_prototype.csv'))};assert len(source)==len(orig)==22
rows=[]
for r in source:
 o=orig[r['variant']];a=r['higher_mvp_insomnia_allele'];b=r['comparison_allele'];assert {a,b}=={o['counted_allele'],o['other_allele']} and float(r['hfrs_beta'])>0 and float(r['mvp_insomnia_or_approx'])>1
 rows.append({'variant':r['variant'],'chr_GRCh37_assumed':o['chr'],'position_GRCh37_assumed':o['position_GRCh37_assumed'],'counted_allele':a,'other_allele':b,'weight':1,'nearest_gene_body_positional':r['nearest_gene_body_positional'],'reference_alleles':r['reference_alleles'],'status':'Exploratory unweighted shared-direction allele count; postselected; no outcome validation or clinical interpretation'})
p=pathlib.Path('outputs/unvalidated_22_shared_direction_score.csv')
with p.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
print('Created',len(rows),'variants. sha256',hashlib.sha256(p.read_bytes()).hexdigest())
