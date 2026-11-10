import requests,csv,json,hashlib
u='https://www.ebi.ac.uk/gwas/rest/api/v2/studies';p={'efo_id':'EFO_0009885','show_child_traits':'true','size':500};r=requests.get(u,params=p,timeout=60);r.raise_for_status();j=r.json();assert j['page']['totalPages']==1
raw='outputs/stage6_frailty_efo_catalog_response.json';open(raw,'w').write(json.dumps(j,indent=2)+'\n');studies=j.get('_embedded',{}).get('studies',[]);rows=[]
for x in studies:rows.append({'accession':x['accession_id'],'trait':x['disease_trait'],'pubmed':x.get('pubmed_id'),'catalog_initial_sample_description':x.get('initial_sample_size'),'catalog_cohorts':';'.join(x.get('cohort',[])),'catalog_full_summary_available':x['full_summary_stats_available'],'summary_url':x.get('full_summary_stats','')})
with open('outputs/stage6_frailty_efo_catalog_inventory.csv','w') as f:
 w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
print('query',r.url,'total',len(rows),'full',sum(r['catalog_full_summary_available'] for r in rows),'response_sha256',hashlib.sha256(open(raw,'rb').read()).hexdigest())
