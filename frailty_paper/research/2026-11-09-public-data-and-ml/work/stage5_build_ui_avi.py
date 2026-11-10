import json,csv,urllib.parse
values=json.load(open('work/stage5_atlas_ui_avi_values.json'));q=[x for x in csv.DictReader(open('outputs/stage5_score_queue.csv')) if 'exploratory_marker' in x['input_classes']];assert len(q)==len(values)==23 and {x['variant_key'] for x in q}==set(values)
r=[]
for x in q:
 v=x['variant_key'];p,top=values[v];r.append({'marker':'|'.join(z for z in x['input_labels'].split('|') if z.startswith('rs')),'variant_GRCh38_ref_alt':v,'avi_phred_displayed_rounded':p,'top_genomewide_snv_percent_displayed':top,'source_url':'https://deepmind.google.com/science/alphagenome/atlas?q='+urllib.parse.quote(v,safe=''),'source':'AlphaGenome Atlas web interface','access_date_utc':'2026-10-27','warning':'display precision only; general variant-impact prediction not phenotype risk; no tissue, gene, mechanism, disease causality or individual effect claim'})
with open('outputs/stage5_23_markers_atlas_avi.csv','w',newline='') as f:w=csv.DictWriter(f,fieldnames=r[0]);w.writeheader();w.writerows(r)
print('rows',len(r),'unique',len({x['variant_GRCh38_ref_alt'] for x in r}))
