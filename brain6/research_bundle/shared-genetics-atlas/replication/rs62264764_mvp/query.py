import collections,csv,gzip,hashlib,json,math,pathlib
p=pathlib.Path(__file__).resolve().parent;root=p.parents[1];prior=root/'replication/mvp_six_leads/receipt.json';mvp=pathlib.Path(json.loads(prior.read_text())['mvp_path']);context=root/'variant_context/selected_variant_gene_context.tsv'
def sha(q):
 h=hashlib.sha256()
 with q.open('rb') as f:
  for x in iter(lambda:f.read(8*1024*1024),b''):h.update(x)
 return h.hexdigest()
source_sha=sha(mvp);assert source_sha==json.loads(prior.read_text())['mvp_sha256']
with context.open(newline='') as f:c=[x for x in csv.DictReader(f,delimiter='\t') if x['snp']=='rs62264764']
assert len(c)==1 and c[0]['source_A1']=='A' and c[0]['grch38_ref_alt']=='chr3:117920728:G>A'
found=[];n=0
with gzip.open(mvp,'rt') as f:
 head=next(f).rstrip('\n').split('\t');ci=head.index('chromosome');pi=head.index('base_pair_location')
 for line in f:
  n+=1;v=line.rstrip('\n').split('\t')
  if v[ci]=='3' and v[pi]=='117920728':found.append(dict(zip(head,v)))
assert n==19703815
matched=[x for x in found if x['rsid']=='rs62264764' and {x['effect_allele'],x['other_allele']}=={'A','G'}];status='NOT_ESTIMABLE';result={'variant':'rs62264764','coordinate':'3:117920728:G>A','aligned_allele':'A','status':status}
if len(matched)==1:
 x=matched[0];odds,pv,freq=map(float,(x['odds_ratio'],x['p_value'],x['effect_allele_frequency']));assert math.isfinite(odds) and odds>0 and math.isfinite(pv) and 0<=pv<=1 and math.isfinite(freq) and 0<=freq<=1
 beta=math.log(odds)*(1 if x['effect_allele']=='A' else -1);discovery=float(c[0]['insomnia_beta_A1']);result.update(status='MATCHED',mvp_p=pv,mvp_logodds_aligned=beta,same_direction_as_original=beta*discovery>0)
elif len(matched)>1:result['status']='AMBIGUOUS'
with (p/'mvp_result.tsv').open('w',newline='') as h:w=csv.DictWriter(h,fieldnames=result.keys(),delimiter='\t');w.writeheader();w.writerow(result)
ld=list(csv.DictReader((p/'pairwise_ld.tsv').open(),delimiter='\t'));assert len(ld)==1
r=ld[0];discordant=round(int(r['haplotypes'])*(float(r['allele1_frequency'])+float(r['allele2_frequency'])-2*float(r['p11'])))
receipt={'mvp_source_sha256':source_sha,'mvp_rows':n,'target_rows_at_position':len(found),'exact_matches':len(matched),'exact_records':matched,'reference_haplotype_discordance':discordant,'reference_haplotypes':int(r['haplotypes']),'note':'Exact allele association and reference-panel matching do not validate causal variant or sample independence.'};(p/'query_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({'result':result,'discordant_reference_haplotypes':discordant},indent=2))
