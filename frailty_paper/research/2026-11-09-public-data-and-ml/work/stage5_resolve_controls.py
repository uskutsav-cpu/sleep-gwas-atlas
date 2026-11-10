import csv,collections,json,subprocess,concurrent.futures,hashlib,pathlib
r=list(csv.DictReader(open('outputs/stage5_matched_control_candidates_grch37.csv')));extra='outputs/stage5_relaxed_control_candidates_grch37.csv'
try:
 for x in csv.DictReader(open(extra)):
  if x['control_rsid']:
   for i,z in enumerate(r):
    if z['candidate']==x['candidate'] and not z['control_rsid']:r[i]=x
except FileNotFoundError:pass
comp=str.maketrans('ACGT','TGCA')
def get(x):
 if not x['control_rsid']:return {**x,'status_GRCh38':'unresolved_no_control'}
 c=x['chr37'];p=x['control_pos37'];rs=x['control_rsid'];u=f'https://grch37.rest.ensembl.org/map/human/GRCh37/{c}:{p}..{p}/GRCh38?content-type=application/json';v='https://rest.ensembl.org/variation/human/'+rs+'?content-type=application/json'
 try:
  a=subprocess.run(['curl','-fsS','--max-time','35',u],capture_output=True,check=True).stdout;maps=json.loads(a).get('mappings',[])
  if len(maps)!=1:raise ValueError('liftover not unique')
  z=maps[0]['mapped'];chrom=str(z['seq_region_name']);pos=int(z['start']);strand=int(z['strand']);assert z['start']==z['end'] and chrom.isdigit()
  alleles=set(x['control_bim_alleles'].split('/'));alleles={q.translate(comp) for q in alleles} if strand==-1 else alleles
  a2=subprocess.run(['curl','-fsS','--max-time','35',v],capture_output=True,check=True).stdout;j=json.loads(a2);matches=[q for q in j.get('mappings',[]) if q['assembly_name']=='GRCh38' and str(q['seq_region_name'])==chrom and q['start']==pos and q['end']==pos and alleles.issubset(set(q['allele_string'].split('/')))];assert len(matches)==1,'no unique rsID plus allele pair at mapped position'
  u3=f'https://rest.ensembl.org/sequence/region/human/{chrom}:{pos}..{pos}:1?content-type=text/plain';a3=subprocess.run(['curl','-fsS','--max-time','35',u3],capture_output=True,check=True).stdout;ref=a3.decode().strip().upper();assert ref in alleles;alt=next(iter(alleles-{ref}));return {**x,'chr_GRCh38':'chr'+chrom,'pos_GRCh38':pos,'ref_GRCh38':ref,'alt_GRCh38':alt,'liftover_strand':strand,'mapping_url':u,'mapping_sha256':hashlib.sha256(a).hexdigest(),'variation_url':v,'variation_sha256':hashlib.sha256(a2).hexdigest(),'reference_url':u3,'reference_sha256':hashlib.sha256(a3).hexdigest(),'status_GRCh38':'unique_liftover_rsid_pair_reference_verified'}
 except Exception as e:return {**x,'status_GRCh38':'needs_review','error':str(e)}
previous=pathlib.Path('outputs/stage5_matched_controls.csv');saved={x['candidate']:x for x in csv.DictReader(previous.open()) if x['status_GRCh38']=='unique_liftover_rsid_pair_reference_verified'} if previous.exists() else {}
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as e:out=list(e.map(lambda x:saved.get(x['candidate']) or get(x),r))
fields=list(r[0])+['chr_GRCh38','pos_GRCh38','ref_GRCh38','alt_GRCh38','liftover_strand','mapping_url','mapping_sha256','variation_url','variation_sha256','reference_url','reference_sha256','status_GRCh38','error']
with open('outputs/stage5_matched_controls.csv','w',newline='') as f:w=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore');w.writeheader();w.writerows(out)
print(dict(collections.Counter(x['status_GRCh38'] for x in out)),flush=True)
