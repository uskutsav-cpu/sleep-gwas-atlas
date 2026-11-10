import csv,json,urllib.parse
score=json.load(open('work/stage5_atlas_ui_credible_avi_values.json'));queue={x['variant_key']:x for x in csv.DictReader(open('outputs/stage5_score_queue.csv'))};pro=list(csv.DictReader(open('outputs/stage5_molecular_credible_set_provenance.csv')));rows=[]
for v,(avi,top) in score.items():
 x=queue[v];src=set(x['input_labels'].split('|'));q=[z for z in pro if z['variant'] in src];assert q
 rows.append({'variant_GRCh38_ref_alt':v,'credible_source_ids':'|'.join(sorted(src)),'rsids':'|'.join(sorted({z['rsid'] for z in q if z['rsid']})),'gene_ids':'|'.join(sorted({z['gene_id'] for z in q if z['gene_id']})),'molecular_pips':'|'.join(sorted({z['molecular_pip'] for z in q if z['molecular_pip']})),'credible_sets':'|'.join(sorted({z['credible_set_id'] for z in q if z['credible_set_id']})),'avi_phred_displayed_rounded':avi,'top_genomewide_snv_percent_displayed':top,'source_url':'https://deepmind.google.com/science/alphagenome/atlas?q='+urllib.parse.quote(v,safe=''),'warning':'molecular PIP and AVI quantify different things; neither is a disease causal probability; website values rounded'})
with open('outputs/stage5_priority_credible_snv_atlas_avi.csv','w',newline='') as f:w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
print('priority alternatives',len(rows))
