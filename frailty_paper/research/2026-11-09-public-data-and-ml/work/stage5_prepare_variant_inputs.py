import csv,collections,pathlib,re,json,urllib.request,hashlib,time,subprocess
files=['outputs/stage3_adarb1_muscle_splice_credible_set.csv','outputs/stage3_skeletal_qtl_credible_set_variants.csv','outputs/stage4_cd4_smad5_credible_set.csv','outputs/stage3_inka1_muscle_credible_set.csv','outputs/stage3_remaining_molecular_credible_sets.csv','outputs/stage3_rs323509_testis_splice_credible_set.csv','outputs/stage3_foxp2_pituitary_susie_credible_set.csv','outputs/stage3_smim32_muscle_credible_set.csv']
records=collections.defaultdict(list)
for f in files:
 for r in csv.DictReader(open(f)):
  v=r['variant'];m=re.fullmatch(r'chr(\d+)_(\d+)_([ACGT]+)_([ACGT]+)',v);assert m,(f,v)
  records[v].append({'source_file':f,'candidate':r.get('candidate',''),'dataset':r.get('dataset',''),'gene_id':r.get('gene_id',''),'molecular_trait_id':r.get('molecular_trait_id',''),'credible_set_id':r.get('cs_id',''),'rsid':r.get('rsid',''),'molecular_pip':r.get('pip',''),'qtl_p':r.get('pvalue','')})
markers=[]
for r in csv.DictReader(open('outputs/shared22_eight_source_evidence_matrix.csv')):
 markers.append({'marker':r['variant'],'chromosome':r['chr_GRCh38'],'position':int(r['position_GRCh38']),'study_counted_allele':r['score_counted_allele'],'study_other_allele':r['other_allele'],'prior_gene':r['nearest_gene_body_positional']})
x=next(x for x in csv.DictReader(open('outputs/integrated_15_candidate_review.csv')) if x['variant']=='rs323509')
markers.append({'marker':'rs323509','chromosome':x['GRCh37_chr'],'position':int(x['GRCh38_pos']),'study_counted_allele':x['counted_allele'],'study_other_allele':x['other_allele'],'prior_gene':x['nearest_gene_body_positional']});assert len(markers)==23
queries={}
for x in markers:queries[(x['chromosome'],x['position'],1)]=None
for v in records:
 c,pos,ref,alt=v.split('_');queries[(c[3:],int(pos),len(ref))]=None
batch=[];keys=list(queries);url='https://rest.ensembl.org/sequence/region/human';failed=[]
for i in range(0,len(keys),25):
 block=keys[i:i+25];regions=[f'{c}:{pos}..{pos+length-1}:1' for c,pos,length in block];payload=json.dumps({'regions':regions}).encode()
 try:
  body=subprocess.run(['curl','--fail','--silent','--show-error','--max-time','45','-X','POST',url,'-H','Content-Type: application/json','-H','Accept: application/json','--data-binary','@-'],input=payload,capture_output=True,check=True).stdout
  data=json.loads(body);mapping={v['query']:v['seq'].upper() for v in data};assert set(mapping)==set(regions)
  for k,q in zip(block,regions):queries[k]=mapping[q]
  batch.append({'n':len(block),'regions':regions,'response_sha256':hashlib.sha256(body).hexdigest(),'status':'ok'})
 except Exception as e:
  failed.append({'regions':regions,'error':str(e)});batch.append({'n':len(block),'regions':regions,'status':'failed','error':str(e)})
 print('batch',i//25+1,'of',(len(keys)+24)//25,'failed',len(failed),flush=True);time.sleep(.15)
rows=[]
for x in markers:
 ref=queries[(x['chromosome'],x['position'],1)];alleles={x['study_counted_allele'],x['study_other_allele']};alt=next(iter(alleles-{ref})) if ref in alleles and len(alleles)==2 else ''
 rows.append({'input_class':'exploratory_marker','input_label':x['marker'],'chr':'chr'+x['chromosome'],'position_GRCh38':x['position'],'reference_allele':ref or '', 'alternate_allele':alt,'status':'reference_verified_study_pair' if alt else 'unresolved_reference_or_study_alleles','study_counted_allele':x['study_counted_allele'],'study_other_allele':x['study_other_allele'],'source_variant_id':'','molecular_records':0,'aliases':x['marker'],'gene_context':x['prior_gene'],'caveat':'GRCh38 marker mapping inherited; sequence check does not independently validate liftover or disease effect alleles'})
for v,es in sorted(records.items()):
 c,pos,ref,alt=v.split('_');observed=queries[(c[3:],int(pos),len(ref))];status='reference_matches_GRCh38' if observed==ref else ('reference_lookup_failed' if observed is None else 'reference_mismatch')
 rows.append({'input_class':'molecular_credible_set_variant','input_label':v,'chr':c,'position_GRCh38':pos,'reference_allele':ref,'alternate_allele':alt,'status':status,'study_counted_allele':'','study_other_allele':'','source_variant_id':v,'molecular_records':len(es),'aliases':'|'.join(sorted({x['rsid'] for x in es if x['rsid']})),'gene_context':'|'.join(sorted({x['gene_id'] for x in es if x['gene_id']})),'caveat':'molecular PIP, not disease PIP; same variant can belong to multiple sets; indel left normalization not yet verified'})
with open('outputs/stage5_variant_input_manifest.csv','w',newline='') as f:w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
links=[{'variant':v,**e} for v,es in records.items() for e in es]
with open('outputs/stage5_molecular_credible_set_provenance.csv','w',newline='') as f:w=csv.DictWriter(f,fieldnames=links[0]);w.writeheader();w.writerows(links)
json.dump({'reference':'Ensembl REST sequence/region human; source response identifies GRCh38; positive strand','url':url,'positions_queried':len(keys),'batches':batch,'failed_batches':failed,'marker_count':len(markers),'distinct_credible_set_variant_count':len(records),'distinct_credible_set_records':len(links),'status_counts':dict(collections.Counter(x['status'] for x in rows))},open('outputs/stage5_variant_reference_check_receipt.json','w'),indent=2)
print('status',dict(collections.Counter(x['status'] for x in rows)),'failed',len(failed),flush=True)
