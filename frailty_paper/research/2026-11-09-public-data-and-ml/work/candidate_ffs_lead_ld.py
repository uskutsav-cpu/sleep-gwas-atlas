import csv,subprocess,json,hashlib
v=list(csv.DictReader(open('outputs/candidate_risk_gene_review_15.csv'))); l=list(csv.DictReader(open('outputs/published_ffs_37_locus_inventory.csv'))); out=[]
for x in v:
 z=sorted((abs(int(x['position_GRCh37'])-int(y['position_reported'])),y['ffs_lead_variant'],y) for y in l if y['chromosome']==x['chr_GRCh37'])
 r={'candidate':x['variant'],'candidate_GRCh37_chr':x['chr_GRCh37'],'candidate_GRCh37_position':x['position_GRCh37'],'nearest_ffs_lead':'','lead_position_reported':'','distance_bp_reported':'','lead_nearest_gene':'','reference_LD_status':'not_queried_over_500kb_or_no_same_chromosome_lead','reference_r2':'','reference_d_prime':'','query_url':'','response_sha256':''}
 if z:
  dist,lead,y=z[0];r.update(nearest_ffs_lead=lead,lead_position_reported=y['position_reported'],distance_bp_reported=dist,lead_nearest_gene=y['nearest_gene_source_annotation_only'])
  if dist<=500000 and lead.startswith('rs'):
   url='https://grch37.rest.ensembl.org/ld/human/pairwise/'+x['variant']+'/'+lead+'?content-type=application/json;population_name=1000GENOMES:phase_3:EUR';r['query_url']=url
   try:
    data=subprocess.check_output(['curl','-fsS','--max-time','25',url]); res=json.loads(data);r['response_sha256']=hashlib.sha256(data).hexdigest()
    if len(res)==1:r.update(reference_LD_status='returned_rsID_based',reference_r2=res[0]['r2'],reference_d_prime=res[0]['d_prime'])
    elif not res:r['reference_LD_status']='empty_result_not_zero_LD'
    else:r['reference_LD_status']='multiple_records_unresolved'
   except Exception as e:r['reference_LD_status']='request_failed'
 out.append(r);print(x['variant'],r['nearest_ffs_lead'],r['reference_LD_status'],r['reference_r2'],flush=True)
with open('outputs/candidate_nearest_ffs_lead_ld.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=out[0]);w.writeheader();w.writerows(out)
