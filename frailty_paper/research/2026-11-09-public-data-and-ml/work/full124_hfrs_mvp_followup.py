import csv,gzip,hashlib,json,pathlib,collections,math
from pyliftover import LiftOver
base=pathlib.Path('/Volumes/Extreme SSD/sleep-frailty-full-autonomous-v5-20261009/frailty_paper/full_autonomous_research_v5/genetics/new_public_validation/raw'); paths={'HFRS':base/'GCST90707291.tsv','MVP':base/'GCST90475826.tsv.gz'}
orig=list(csv.DictReader(open('outputs/published_insomnia_fi_124_cross_trait_candidates.csv'))); assert len(orig)==len({x['variant_id'] for x in orig})==124
lo=LiftOver('work/hg19ToHg38.over.chain.gz');comp=str.maketrans('ACGT','TGCA');targets={}
for r in orig:
 chrom='chr'+r['chromosome_reported']; pos=int(r['position_reported']); hits=lo.convert_coordinate(chrom,pos-1) or []
 hits=[x for x in hits if x[0]=='chr'+r['chromosome_reported']];r['lift_hits']=hits
 if len(hits)==1:
  c,b,strand,_=hits[0];a=r['effect_allele_reported'];o=r['other_allele_reported'];r.update(mapped_chr=c[3:],mapped_pos=str(b+1),mapped_strand=strand,allele_38=a if strand=='+' else a.translate(comp),other_38=o if strand=='+' else o.translate(comp))
 targets[r['variant_id']]=r
found={};receipts={}
for source,p in paths.items():
 acc=collections.defaultdict(list); n=0; md5=hashlib.md5()
 if source=='HFRS':
  with p.open('rb') as f:
   header=None
   for line in f:
    md5.update(line)
    if header is None:header=line.decode().rstrip('\r\n').split('\t'); idx=header.index('rs_id'); continue
    n+=1;v=line.decode().rstrip('\r\n').split('\t');rs=v[idx]
    if rs in targets:acc[rs].append(dict(zip(header,v)))
  assert md5.hexdigest()=='38d9434f50a4493def9d2c9d1cc1b4b3'
  digest=md5.hexdigest()
 else:
  sha=hashlib.sha256()
  with p.open('rb') as f:
   for block in iter(lambda:f.read(8*1024*1024),b''):sha.update(block)
  assert sha.hexdigest()=='46b62349ba0d98b22bf302be973bb54a80707cfb02ce953426f8224910761835'
  digest=sha.hexdigest()
  with gzip.open(p,'rt') as f:
   rr=csv.DictReader(f,delimiter='\t')
   for row in rr:
    n+=1;rs=row['rsid']
    if rs in targets:acc[rs].append(row)
 found[source]=acc; receipts[source]={'n_rows':n,'source_hash':digest,'id_found':len(acc)};print(source,receipts[source],flush=True)
results=[]
for r in orig:
 v=r['variant_id'];row={'variant':v,'song_cpassoc_rank':r['rank_by_cpassoc_p'],'song_chr':r['chromosome_reported'],'song_position_reported':r['position_reported'],'song_build_status':'not verified; GRCh37 assumed for mapping','MHC_flag':r['approx_chr6_MHC_flag'],'GRCh38_chr':r.get('mapped_chr',''),'GRCh38_pos':r.get('mapped_pos',''),'mapping_strand':r.get('mapped_strand',''),'song_allele_assumed_38':r.get('allele_38',''),'song_other_assumed_38':r.get('other_38','')}
 for src in ('HFRS','MVP'):
  hits=found[src].get(v,[]);poshits=[x for x in hits if x['chromosome']==r.get('mapped_chr') and x['base_pair_location']==r.get('mapped_pos')];allelehits=[x for x in poshits if {x['effect_allele'],x['other_allele']}=={r.get('allele_38'),r.get('other_38')}];row[src+'_id_rows']=len(hits);row[src+'_position_rows']=len(poshits);row[src+'_allele_rows']=len(allelehits)
  if len(allelehits)==1:
   x=allelehits[0];row.update({src+'_effect_allele':x['effect_allele'],src+'_other_allele':x['other_allele'],src+'_p':x['p_value'],src+'_effect':x['beta'] if src=='HFRS' else x['odds_ratio'],src+'_se':x['standard_error']})
  else:row.update({src+'_effect_allele':'',src+'_other_allele':'',src+'_p':'',src+'_effect':'',src+'_se':''})
 ok=(len(r['lift_hits'])==1 and r.get('mapped_strand')=='+' and len(r['allele_38'])==len(r['other_38'])==1 and set([r['allele_38'],r['other_38']])<=set('ACGT') and ''.join(sorted(r['allele_38']+r['other_38'])) not in ('AT','CG') and row['HFRS_allele_rows']==row['MVP_allele_rows']==1)
 row['eligible_for_joint_screen']=ok
 if ok:
  try:
   ph,pm=float(row['HFRS_p']),float(row['MVP_p']);bh=float(row['HFRS_effect']);or_=float(row['MVP_effect']);assert all(map(math.isfinite,(ph,pm,bh,or_))) and 0<=ph<=1 and 0<=pm<=1 and or_>0
   # align MVP to HFRS counted allele; OR rounded and identical to 1 has no usable sign
   m=math.log(or_)*(1 if row['HFRS_effect_allele']==row['MVP_effect_allele'] else -1);row.update(joint_p_max=max(ph,pm),joint_bonferroni_124=min(1,124*max(ph,pm)),same_effect_direction=bh*m>0 if bh*m else 'zero_or_rounded')
  except (ValueError,AssertionError):row['eligible_for_joint_screen']=False
 results.append(row)
keys=list(dict.fromkeys(k for row in results for k in row));dest=pathlib.Path('outputs/song124_mvp_hfrs_joint_screen.csv')
with dest.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(results)
eligible=[r for r in results if r['eligible_for_joint_screen']];sig=[r for r in eligible if r['joint_bonferroni_124']<.05]
receipt={'inputs':{k:str(v) for k,v in paths.items()},'sources':receipts,'source_candidate_sha256':hashlib.sha256(pathlib.Path('outputs/published_insomnia_fi_124_cross_trait_candidates.csv').read_bytes()).hexdigest(),'liftover_chain_sha256':hashlib.sha256(pathlib.Path('work/hg19ToHg38.over.chain.gz').read_bytes()).hexdigest(),'output_sha256':hashlib.sha256(dest.read_bytes()).hexdigest(),'n_candidates':len(orig),'n_eligible':len(eligible),'n_bonferroni124_screen_below_005':len(sig),'n_nominal':sum(r['joint_p_max']<.05 for r in eligible),'threshold_hits':[{'variant':r['variant'],'joint_p':r['joint_p_max'],'same_direction':r['same_effect_direction'],'MHC_flag':r['MHC_flag']} for r in sorted(sig,key=lambda x:x['joint_p_max'])]}
pathlib.Path('outputs/song124_mvp_hfrs_joint_screen.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2),flush=True)
