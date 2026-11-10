import urllib.request,urllib.parse,json,csv,hashlib,time,subprocess
pairs=[('prototype','rs113851554','rs78206187',271670),('prototype','rs11152363','rs1942262',183871),('prototype','rs35881094','rs17550821',55602),('published_ffs_insomnia_region','rs7703746','rs1592757',17339),('song_chr3_candidates','rs34290943','rs10865954',50329)]
urlbase='https://grch37.rest.ensembl.org/ld/human/pairwise/'
rows=[]
for group,a,b,dist in pairs:
 url=urlbase+a+'/'+b+'?content-type=application/json;population_name=1000GENOMES:phase_3:EUR'
 out={'pair_type':group,'variant_1':a,'variant_2':b,'reported_distance_bp':dist,'population':'1000GENOMES:phase_3:EUR','source_url':url,'status':'','r2':'','d_prime':'','interpretation':'Pairwise LD reference only; no causal colocalization, replication or clinical score validation'}
 try:
  body=subprocess.run(["curl","--fail","--silent","--show-error","--max-time","25",url],capture_output=True,check=True).stdout;data=json.loads(body)
  if isinstance(data,list) and len(data)==1:
   out.update(status='returned',r2=data[0]['r2'],d_prime=data[0]['d_prime'],response_sha256=hashlib.sha256(body).hexdigest())
  elif isinstance(data,list) and len(data)==0:out.update(status='no_result',response_sha256=hashlib.sha256(body).hexdigest())
  else:out.update(status='unexpected_response',response_sha256=hashlib.sha256(body).hexdigest())
 except Exception as e:out.update(status='request_failed',response_sha256='',error=str(e))
 rows.append(out); print(a,b,out['status'],out['r2'],flush=True);time.sleep(.15)
keys=['pair_type','variant_1','variant_2','reported_distance_bp','population','status','r2','d_prime','response_sha256','source_url','interpretation']
with open('outputs/selected_variant_pair_ld_reference.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=keys,extrasaction='ignore');w.writeheader();w.writerows(rows)
