import csv,gzip,collections,json,pathlib,hashlib,datetime,math
from pyliftover import LiftOver
p=pathlib.Path('outputs');q=list(csv.DictReader(open(p/'stage5_score_queue.csv')));f=LiftOver('work/hg38ToHg19.over.chain.gz');back=LiftOver('work/hg19ToHg38.over.chain.gz');mapped={};audit=[]
for x in q:
 k=x['variant_key'];ref=x['reference_allele'];alt=x['alternate_allele'];entry={'variant_GRCh38':k,'status':'','chr37':'','pos37':'','source_ref':ref,'source_alt':alt}
 if len(ref)!=1 or len(alt)!=1:entry['status']='indel_not_attempted'
 elif {ref,alt} in ({'A','T'},{'C','G'}):entry['status']='palindromic_snv_excluded'
 else:
  candidates=f.convert_coordinate(x['chr'],int(x['position_GRCh38'])-1)
  if len(candidates)!=1 or candidates[0][2]!='+':entry['status']='no_unique_positive_forward_mapping'
  else:
   chrom,pos,strand,_=candidates[0];reverse=back.convert_coordinate(chrom,pos)
   if len(reverse)!=1 or reverse[0][:3]!=(x['chr'],int(x['position_GRCh38'])-1,'+'):entry['status']='reciprocal_mapping_failed'
   else:
    entry.update(status='unique_positive_reciprocal_mapping',chr37=chrom.removeprefix('chr'),pos37=str(pos+1));mapped[(chrom.removeprefix('chr'),str(pos+1),frozenset((ref,alt)))]=x
 audit.append(entry)
with open(p/'stage5_queue_snv_reciprocal_liftover.csv','w',newline='') as out:
 w=csv.DictWriter(out,fieldnames=list(audit[0]));w.writeheader();w.writerows(audit)
traits=['Atkins_FI','Jansen_insomnia','Mak_FI','Mak_FP','Ye_FFS'];matches=collections.defaultdict(list);source_info={}
for t in traits:
 path=pathlib.Path('work/stage3_region_cache')/(t+'.regional.tsv.gz');source_info[t]=hashlib.sha256(path.read_bytes()).hexdigest()
 with gzip.open(path,'rt') as inp:
  for row in csv.DictReader(inp,delimiter='\t'):
   chrom=row.get('CHR',row.get('CHROM'));pos=row.get('BP',row.get('POS'));a1=row.get('A1','');a2=row.get('A2','');x=mapped.get((chrom,pos,frozenset((a1,a2))))
   if x is None or a1==a2:continue
   beta=float(row['BETA']);alt=x['alternate_allele'];af=row.get('FRQ',row.get('AF1',''));alt_af='' if not af else str(float(af) if a1==alt else 1-float(af));matches[(x['variant_key'],t)].append({'beta':beta if a1==alt else -beta,'se':row['SE'],'p':row['P'],'alt_af':alt_af,'a1':a1,'a2':a2,'source_rsid':row.get('SNP',''),'source_chr37':chrom,'source_pos37':pos,'n':row.get('N','')})
rows=[];statusmap={x['variant_GRCh38']:x['status'] for x in audit}
for x in q:
 for t in traits:
  key=x['variant_key'];hits=matches.get((key,t),[]);status=statusmap[key]
  if status=='unique_positive_reciprocal_mapping':status='unique_allele_match' if len(hits)==1 else 'no_exact_allele_match_in_cache' if not hits else 'multiple_matching_rows_excluded'
  h=hits[0] if status=='unique_allele_match' else {}
  rows.append({'variant_GRCh38':key,'trait_source':t,'status':status,'chr37':h.get('source_chr37',''),'pos37':h.get('source_pos37',''),'source_variant_id':h.get('source_rsid',''),'source_effect_allele_A1':h.get('a1',''),'source_other_allele_A2':h.get('a2',''),'gwas_beta_alt':h.get('beta',''),'gwas_se':h.get('se',''),'gwas_p':h.get('p',''),'gwas_af_alt':h.get('alt_af',''),'gwas_n':h.get('n',''),'source_cache_sha256':source_info[t],'caveat':'GRCh37 source build assumed; positive reciprocal chain and allele-pair match, not independent reference sequence confirmation; Ye A1 effect direction inferred; no causal inference'})
with open(p/'stage5_grch37_gwas_exact_queue.csv','w',newline='') as out:w=csv.DictWriter(out,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
existing=list(csv.DictReader(open(p/'stage5_credible_variants_gwas_aligned_evidence.csv')));table={(r['variant_GRCh38'],r['trait_source']):r for r in rows};checked=0;discrepancies=[]
for r in existing:
 if r['gwas_trait_source'] not in traits:continue
 x=table[(r['variant_GRCh38'],r['gwas_trait_source'])]
 if x['status']!='unique_allele_match':continue
 checked+=1
 if abs(float(x['gwas_beta_alt'])-float(r['gwas_beta_alt']))>1e-8 or abs(float(x['gwas_p'])-float(r['gwas_p']))>1e-8:discrepancies.append((r['variant_GRCh38'],r['gwas_trait_source']))
summary={'date_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'queue':len(q),'mapping_status':dict(collections.Counter(x['status'] for x in audit)),'match_status_by_trait':{t:dict(collections.Counter(x['status'] for x in rows if x['trait_source']==t)) for t in traits},'prior_aligned_rows_checked':checked,'prior_discrepancies':discrepancies,'chain38to37_sha256':hashlib.sha256(pathlib.Path('work/hg38ToHg19.over.chain.gz').read_bytes()).hexdigest(),'chain37to38_sha256':hashlib.sha256(pathlib.Path('work/hg19ToHg38.over.chain.gz').read_bytes()).hexdigest(),'pyliftover_version':'0.4.1','caveat':'An exact allele-pair match with reciprocal chain does not independently verify strand or source build, or validate association.'}
(p/'stage5_grch37_gwas_exact_queue_audit.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2));assert not discrepancies
