import csv,json,hashlib,urllib.request,ssl,time
ctx=ssl._create_unverified_context()
rows=list(csv.DictReader(open('outputs/prototype_variant_identity_reference.csv')))
out=[]
for r in rows:
 url='https://grch37.rest.ensembl.org/vep/human/id/'+r['variant']+'?content-type=application/json'
 try:
  data=urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'research-audit/1.0'}),context=ctx,timeout=30).read(); hits=json.loads(data)
  hit=[x for x in hits if str(x['seq_region_name'])==r['chr_expected'] and str(x['start'])==r['pos_expected']]
  if len(hit)!=1: raise ValueError('coordinate records '+str(len(hit)))
  h=hit[0]; alleles=h['allele_string'].split('/'); ref=alleles[0]; pair={r['A1'],r['A2']}
  if ref not in pair or len(pair)!=2: raise ValueError('study pair is not reference-alternate')
  alt=next(iter(pair-{ref}))
  if alt not in alleles: raise ValueError('alternate absent')
  tc=[x for x in h.get('transcript_consequences',[]) if x.get('variant_allele')==alt]
  terms=sorted({term for x in tc for term in x.get('consequence_terms',[])})
  genes=sorted({x.get('gene_symbol') or x.get('gene_id') for x in tc if x.get('gene_symbol') or x.get('gene_id')})
  regulatory=sorted({term for x in h.get('regulatory_feature_consequences',[]) if x.get('variant_allele')==alt for term in x.get('consequence_terms',[])})
  out.append({'variant':r['variant'],'chr_GRCh37':r['chr_expected'],'position_GRCh37':r['pos_expected'],'counted_allele':r['A1'],'other_allele':r['A2'],'reference_allele':ref,'contrast_alternate_allele':alt,'vep_all_alleles':h['allele_string'],'transcript_terms':';'.join(terms),'transcript_genes':';'.join(genes),'regulatory_terms':';'.join(regulatory),'transcript_count':len(tc),'status':'Allele-specific VEP transcript annotations; no disease mechanism established','source_url':url,'response_sha256':hashlib.sha256(data).hexdigest()})
  print(r['variant'],alt,','.join(terms),flush=True)
 except Exception as e: print(r['variant'],'ERROR',str(e),flush=True)
with open('outputs/candidate_vep_allele_consequences.csv','w',newline='') as f:
 if out:
  w=csv.DictWriter(f,fieldnames=list(out[0]));w.writeheader();w.writerows(out)
print('complete',len(out),'of',len(rows))
