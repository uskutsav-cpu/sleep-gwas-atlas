import requests,csv,hashlib,json,math
u='https://www.nhlbi.nih.gov/grasp/downloads/FullResults/2017/2017_Ben_Avrahim/Gait-Aging-All.txt';c=list(csv.DictReader(open('outputs/expanded_22_candidate_risk_gene_table.csv')));ids={x['variant'] for x in c};frequency={x['variant']:x for x in csv.DictReader(open('outputs/22_score_reference_group_allele_frequencies.csv'))};found={k:[] for k in ids};h=hashlib.sha256();n=0
with requests.get(u,stream=True,timeout=(30,180)) as r:
 r.raise_for_status();it=r.iter_lines();head=next(it);h.update(head+b'\n');columns=head.decode().split();assert columns[0]=='all' and columns[2:4]==['Allele1','Allele2'];
 for line in it:
  if not line:continue
  h.update(line+b'\n');parts=line.decode().split();assert len(parts)==len(columns),(n,len(parts));n+=1
  if parts[0] in found:found[parts[0]].append(dict(zip(columns,parts)))
out=[]
for x in c:
 rows=found[x['variant']];a=x['counted_allele'].upper();b=x['other_allele'].upper();m=[z for z in rows if {z['Allele1'].upper(),z['Allele2'].upper()}=={a,b}];d={'variant':x['variant'],'candidate_allele':a,'candidate_other':b,'source_rsid_rows':len(rows),'source_matching_allele_pairs':len(m),'status':'not_found' if not rows else 'allele_mismatch' if not m else 'duplicate' if len(m)>1 else 'matched','palindromic':{a,b} in [{'A','T'},{'C','G'}]}
 if len(m)==1:
  z=m[0];sgn=1 if z['Allele1'].upper()==a else -1;beta=sgn*float(z['Effect']);af=float(z['Freq1']);af=af if sgn==1 else 1-af;f=float(frequency[x['variant']]['reference_EUR_frequency']);d.update({'source_allele1':z['Allele1'],'source_allele2':z['Allele2'],'beta_aligned':beta,'se':z['StdErr'],'p':z['P-value'],'N_weight':z['Weight'],'candidate_af':af,'EUR_reference_af':f,'absolute_af_difference':abs(af-f),'direction':z['Direction'],'heterogeneity_p':z['HetPVal']})
 out.append(d)
fields=list(dict.fromkeys(k for d in out for k in d));
with open('outputs/stage6_charge_gait_full_candidate_lookup.csv','w') as f:
 w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(out)
receipt={'source_url':u,'decoded_rows':n,'decoded_response_sha256_lf_normalized':h.hexdigest(),'response_etag':r.headers.get('ETag'),'last_modified':r.headers.get('Last-Modified'),'candidate_count':len(c),'status_counts':{s:sum(z['status']==s for z in out) for s in ['matched','not_found','allele_mismatch','duplicate']},'warning':'rsID and allele pair matching without study-specific genomic coordinates does not resolve strand or multiallelic ambiguity'}
open('outputs/stage6_charge_gait_full_lookup_receipt.json','w').write(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
