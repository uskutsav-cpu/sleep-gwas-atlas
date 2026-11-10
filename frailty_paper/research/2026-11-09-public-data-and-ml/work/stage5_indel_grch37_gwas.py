import csv,json,collections,concurrent.futures,subprocess,hashlib,pathlib,gzip,datetime
from pyliftover import LiftOver
out=pathlib.Path('outputs');q=[x for x in csv.DictReader(open(out/'stage5_score_queue.csv')) if len(x['reference_allele'])!=len(x['alternate_allele'])];f=LiftOver('work/hg38ToHg19.over.chain.gz');b=LiftOver('work/hg19ToHg38.over.chain.gz');cache=pathlib.Path('work/stage5_indel_hg19_ref');cache.mkdir(exist_ok=True)
def run(x):
 k=x['variant_key'];chrom=x['chr'];pos=int(x['position_GRCh38']);ref=x['reference_allele'];alt=x['alternate_allele'];m=[]
 for offset in range(len(ref)+1):
  y=f.convert_coordinate(chrom,pos-1+offset)
  if len(y)!=1 or y[0][2]!='+':return {'variant_GRCh38':k,'status':'mapping_not_unique_positive'}
  z=y[0];back=b.convert_coordinate(z[0],z[1]);
  if len(back)!=1 or back[0][:3]!=(chrom,pos-1+offset,'+'):return {'variant_GRCh38':k,'status':'reciprocal_failure'}
  m.append(z)
 if any(y[0]!=m[0][0] or y[1]!=m[0][1]+i for i,y in enumerate(m)):return {'variant_GRCh38':k,'status':'reference_span_not_contiguous'}
 c=m[0][0];pos37=m[0][1]+1;start=max(0,pos37-101);end=pos37+len(ref)+100;url=f'https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom={c};start={start};end={end}';result=subprocess.run(['curl','-fsS','--max-time','30',url],capture_output=True,check=True);j=json.loads(result.stdout);seq=j['dna'].upper();assert j['start']==start and len(seq)==end-start
 id=hashlib.sha256(k.encode()).hexdigest()[:20];(cache/(id+'.json')).write_text(json.dumps({'variant':k,'url':url,'response':j}))
 if seq[pos37-1-start:pos37-1-start+len(ref)]!=ref:return {'variant_GRCh38':k,'status':'hg19_reference_mismatch','hg19_position':pos37,'sequence_url':url}
 original=seq[:pos37-1-start]+alt+seq[pos37-1-start+len(ref):];n=pos37;r=ref;a=alt
 while r[-1]==a[-1]:
  if min(len(r),len(a))==1:
   if n-1<=start:return {'variant_GRCh38':k,'status':'normalization_window_exhausted'}
   prev=seq[n-2-start];r=prev+r;a=prev+a;n-=1
  r=r[:-1];a=a[:-1]
 while len(r)>1 and len(a)>1 and r[0]==a[0]:r=r[1:];a=a[1:];n+=1
 assert seq[:n-1-start]+a+seq[n-1-start+len(r):]==original
 return {'variant_GRCh38':k,'status':'mapped_reference_checked_normalized','hg19_chrom':c.removeprefix('chr'),'hg19_position':n,'hg19_ref':r,'hg19_alt':a,'shift_37':n-pos37,'sequence_url':url}
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool: result=list(pool.map(run,q))
fields=['variant_GRCh38','status','hg19_chrom','hg19_position','hg19_ref','hg19_alt','shift_37','sequence_url']
with open(out/'stage5_indel_hg19_mapping.csv','w',newline='') as z:w=csv.DictWriter(z,fieldnames=fields);w.writeheader();w.writerows(result)
lookup={(x['hg19_chrom'],str(x['hg19_position']),frozenset((x['hg19_ref'],x['hg19_alt']))):x for x in result if x['status']=='mapped_reference_checked_normalized'}
traits=['Atkins_FI','Jansen_insomnia','Mak_FI','Mak_FP','Ye_FFS'];matches=collections.defaultdict(list)
for t in traits:
 path=pathlib.Path('work/stage3_region_cache')/(t+'.regional.tsv.gz')
 with gzip.open(path,'rt') as z:
  for x in csv.DictReader(z,delimiter='\t'):
   key=(x.get('CHR',x.get('CHROM')),x.get('BP',x.get('POS')),frozenset((x['A1'],x['A2'])));y=lookup.get(key)
   if not y or x['A1']==x['A2']:continue
   v=float(x['BETA']);matches[(y['variant_GRCh38'],t)].append({'source_id':x.get('SNP',''),'beta_alt':v if x['A1']==y['hg19_alt'] else -v,'se':x['SE'],'p':x['P'],'source_a1':x['A1'],'source_a2':x['A2']})
rows=[]
for x in result:
 for t in traits:
  z=matches.get((x['variant_GRCh38'],t),[]);status=x['status'] if x['status']!='mapped_reference_checked_normalized' else 'unique_exact_allele_match' if len(z)==1 else 'no_exact_representation_match' if len(z)==0 else 'multiple_rows_excluded';r=z[0] if status=='unique_exact_allele_match' else {}
  rows.append({'variant_GRCh38':x['variant_GRCh38'],'trait_source':t,'match_status':status,**{v:r.get(v,'') for v in ['source_id','beta_alt','se','p','source_a1','source_a2']},'caveat':'GRCh37 source build assumed; only normalized exact source representation included; source A1 effect coding and Ye direction require confirmation; no causal inference'})
with open(out/'stage5_indel_grch37_gwas_matches.csv','w',newline='') as z:w=csv.DictWriter(z,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
audit={'date_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'indels':len(q),'mapping_status':dict(collections.Counter(x['status'] for x in result)),'source_match_counts':{t:dict(collections.Counter(x['match_status'] for x in rows if x['trait_source']==t)) for t in traits},'positions_left_shifted_in_hg19':sum(x.get('shift_37',0)!=0 for x in result),'note':'Nonmatched indels may be absent or use an alternative representation; does not imply null effects.'};(out/'stage5_indel_grch37_gwas_audit.json').write_text(json.dumps(audit,indent=2)+'\n');print(json.dumps(audit,indent=2))
