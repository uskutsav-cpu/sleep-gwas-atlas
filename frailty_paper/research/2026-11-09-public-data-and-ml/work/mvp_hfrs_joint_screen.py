import csv, pathlib, hashlib
src=pathlib.Path('outputs/integrated_15_candidate_review.csv'); rows=list(csv.DictReader(src.open())); assert len(rows)==15 and len({r['variant'] for r in rows})==15
result=[]
for r in rows:
 p1=float(r['MVP_insomnia_p']); p2=float(r['FinnGen_HFRS_p']); b1=float(r['MVP_insomnia_beta_on_counted_allele']); b2=float(r['FinnGen_HFRS_beta_on_counted_allele']); assert 0<=p1<=1 and 0<=p2<=1 and b1 and b2
 j=max(p1,p2)
 result.append({'variant':r['variant'],'counted_allele':r['counted_allele'],'other_allele':r['other_allele'],'positional_gene_label':r['nearest_gene_body_positional'],'mvp_insomnia_beta':b1,'mvp_insomnia_p':p1,'finn_hfrs_beta':b2,'finn_hfrs_se':r['FinnGen_HFRS_se'],'finn_hfrs_p':p2,'same_sign':b1*b2>0,'joint_p_max':j,'joint_p_times_15_capped_at_1':min(1,15*j),'flags':r['flags']})
result.sort(key=lambda x:x['joint_p_max'])
out=pathlib.Path('outputs/mvp_insomnia_hfrs_joint_candidate_screen.csv')
with out.open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(result[0])); w.writeheader(); w.writerows(result)
print('input_sha256',hashlib.sha256(src.read_bytes()).hexdigest()); print('output_sha256',hashlib.sha256(out.read_bytes()).hexdigest())
for r in result: print(r['variant'],r['joint_p_max'],r['joint_p_times_15_capped_at_1'],r['same_sign'])
