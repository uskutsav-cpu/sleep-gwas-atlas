import csv,json,hashlib,concurrent.futures,subprocess,time
with open('outputs/unvalidated_15_variant_score_prototype.csv',newline='') as f: rows=list(csv.DictReader(f))
def one(c):
 chrom=c['chr'];p=int(c['position_GRCh37_inferred']);lo=max(1,p-100000);hi=p+100000
 url=f'https://grch37.rest.ensembl.org/overlap/region/human/{chrom}:{lo}-{hi}?feature=gene;content-type=application/json'
 for tries in range(3):
  v=subprocess.run(['curl','--fail','--silent','--show-error','--max-time','25',url],capture_output=True)
  if v.returncode==0:
   try:
    genes=json.loads(v.stdout)
    if not isinstance(genes,list):raise ValueError('not list')
    viable=[]
    for x in genes:
     a,b=x.get('start'),x.get('end')
     if not isinstance(a,int) or not isinstance(b,int):continue
     dist=max(a-p,0,p-b)
     viable.append((dist,x.get('id',''),x))
    if not viable:return c,url,hashlib.sha256(v.stdout).hexdigest(),[],''
    m=min(t[0] for t in viable)
    return c,url,hashlib.sha256(v.stdout).hexdigest(),[(d,g,x) for d,g,x in viable if d==m],''
   except Exception as ex: err=str(ex)
  else:err=v.stderr.decode(errors='replace')[:300]
  time.sleep(1+tries)
 return c,url,'',[],err
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex: results=list(ex.map(one,rows))
fields=['variant','chr_GRCh37','position_GRCh37','counted_allele','reference_multiallelic','lookup_status','nearest_gene_ids','nearest_gene_names','nearest_gene_biotypes','distance_to_gene_body_bp','nearest_gene_relation','search_window_bp','ensembl_response_sha256','source_url','interpretation']
with open('outputs/candidate_nearest_gene_annotations.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
 for c,url,h,found,err in results:
  x=dict(variant=c['variant'],chr_GRCh37=c['chr'],position_GRCh37=c['position_GRCh37_inferred'],counted_allele=c['counted_allele'],reference_multiallelic=c['GRCh37_reference_multiallelic'],lookup_status=('error: '+err if err else ('found' if found else 'no_gene_in_search_window')),search_window_bp=100000,ensembl_response_sha256=h,source_url=url,interpretation='Nearest gene body only; does not establish causal gene, function, expression or shared association')
  if found:
   x.update(nearest_gene_ids=';'.join(z[1] for z in found),nearest_gene_names=';'.join(z[2].get('external_name','') for z in found),nearest_gene_biotypes=';'.join(z[2].get('biotype','') for z in found),distance_to_gene_body_bp=found[0][0],nearest_gene_relation='overlaps_gene_body' if found[0][0]==0 else 'nearest_gene_body_within_window')
  w.writerow(x)
print(json.dumps({'n':len(results),'ok':sum(bool(r[3]) for r in results),'empty':sum(not r[3] and not r[4] for r in results),'errors':[r[0]['variant']+': '+r[4] for r in results if r[4]]}))
