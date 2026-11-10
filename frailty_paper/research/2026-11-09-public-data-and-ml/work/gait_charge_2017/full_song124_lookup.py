import csv,hashlib,json,requests
source='https://www.nhlbi.nih.gov/grasp/downloads/FullResults/2017/2017_Ben_Avrahim/Gait-Aging-All.txt';c=list(csv.DictReader(open('outputs/published_insomnia_fi_124_cross_trait_candidates.csv')));ids={x['variant_id'] for x in c};assert len(c)==len(ids)==124;h=hashlib.sha256();hits={k:[] for k in ids};n=0
with requests.get(source,stream=True,timeout=(30,180)) as r:
 r.raise_for_status();lines=r.iter_lines();head=next(lines);h.update(head+b'\n');cols=head.decode().split()
 for line in lines:
  if not line:continue
  h.update(line+b'\n');p=line.decode().split();assert len(p)==len(cols);n+=1
  if p[0] in ids:hits[p[0]].append(dict(zip(cols,p)))
rows=[]
for x in c:
 hts=hits[x['variant_id']];d={'variant':x['variant_id'],'source_rows':len(hts),'status':'unique_rsid' if len(hts)==1 else 'not_found' if not hts else 'duplicate_rsid','paper_effect_allele_not_validated':x['effect_allele_reported'],'paper_other_allele_not_validated':x['other_allele_reported']}
 if len(hts)==1:
  z=hts[0];d.update({'source_allele1':z['Allele1'],'source_allele2':z['Allele2'],'source_p':z['P-value'],'source_N':z['Weight'],'source_effect':z['Effect'],'source_se':z['StdErr'],'source_heterogeneity_p':z['HetPVal']})
 rows.append(d)
with open('outputs/stage6_charge_gait_song124_lookup.csv','w') as f:
 w=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for d in rows for k in d)));w.writeheader();w.writerows(rows)
receipt={'source_url':source,'response_sha256_lf_normalized':h.hexdigest(),'source_rows':n,'candidate_count':len(c),'statuses':{z:sum(d['status']==z for d in rows) for z in ['unique_rsid','not_found','duplicate_rsid']},'note':'rsID matching only; no build, direction, source-allele or multiallelic equivalence established'};open('outputs/stage6_charge_gait_song124_receipt.json','w').write(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
