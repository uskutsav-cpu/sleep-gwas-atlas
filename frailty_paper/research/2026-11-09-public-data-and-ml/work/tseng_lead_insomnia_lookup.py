import csv,pathlib,json,collections,hashlib,math
from pyliftover import LiftOver
lead=list(csv.reader(open('work/tseng_lead_positions.tsv'),delimiter='\t'));pv=dict(csv.reader(open('work/tseng_lead_pvalues.tsv'),delimiter='\t'));lo=LiftOver('work/hg19ToHg38.over.chain.gz')
cols={'MVP':'chromosome base_pair_location effect_allele other_allele odds_ratio standard_error effect_allele_frequency p_value rsid ci_upper ci_lower alt n case_af num_cases control_af num_controls r2 q_pval i2 direction'.split(),'FinnGen':'chrom pos ref alt rsids nearest_genes pval mlogp beta sebeta af_alt af_alt_cases af_alt_controls'.split()}
raw={}
for tag in cols:
 idx=collections.defaultdict(list)
 with open('work/tseng_leads_'+tag.lower()+'_raw.tsv') as f:
  for line in f:
   fields=line.rstrip('\r\n').split('\t'); assert len(fields)==len(cols[tag]);row=dict(zip(cols[tag],fields)); rs=row['rsid'] if tag=='MVP' else row['rsids']
   for t in rs.split(','):idx[t].append(row)
 raw[tag]=idx
out=[]
for rs,ch,pos,a,b in lead:
 rec={'variant':rs,'tseng_chr_assumed_GRCh37':ch,'tseng_pos_assumed_GRCh37':pos,'tseng_a1_reported':a,'tseng_a2_reported':b,'tseng_FI_p':pv[rs]}
 is_snp=len(a)==len(b)==1 and set(a+b)<=set('ACGT');nonpal=is_snp and ''.join(sorted(a+b)) not in ('AT','CG');maps=lo.convert_coordinate('chr'+ch,int(pos)-1) or []; good=[m for m in maps if m[0]=='chr'+ch and m[2]=='+'];rec.update(nonpal_snp=nonpal,liftover_count=len(maps))
 for tag in cols:
  allhits=raw[tag].get(rs,[]);rec[tag+'_rsid_rows']=len(allhits);hits=[]
  if nonpal and len(maps)==len(good)==1:
   target=str(good[0][1]+1)
   for h in allhits:
    c=h['chromosome'] if tag=='MVP' else h['chrom'];p=h['base_pair_location'] if tag=='MVP' else h['pos'];aa=(h['effect_allele'],h['other_allele']) if tag=='MVP' else (h['alt'],h['ref'])
    if c==ch and p==target and set(aa)=={a,b}:hits.append(h)
  rec[tag+'_exact_pair_hits']=len(hits)
  if len(hits)==1:
   h=hits[0]
   if tag=='MVP':
    sign=1 if h['effect_allele']==a else -1; est=float(h['odds_ratio']);beta=sign*math.log(est);prob=float(h['p_value']);qual=h['r2'];freq=h['effect_allele_frequency']
   else:
    sign=1 if h['alt']==a else -1;beta=sign*float(h['beta']);prob=float(h['pval']);qual='';freq=h['af_alt']
   rec.update({tag+'_p':prob,tag+'_beta_or_logOR_aligned_a1':beta,tag+'_source_frequency':freq,tag+'_source_quality_r2':qual})
 if rec['MVP_exact_pair_hits']==rec['FinnGen_exact_pair_hits']==1:
  rec['both_nominal']=rec['MVP_p']<.05 and rec['FinnGen_p']<.05;rec['sleep_beta_agree']=rec['MVP_beta_or_logOR_aligned_a1']*rec['FinnGen_beta_or_logOR_aligned_a1']>0
 out.append(rec)
p=pathlib.Path('outputs/tseng_lead_insomnia_followup.csv');fields=list(dict.fromkeys(k for x in out for k in x))
with p.open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(out)
good=[r for r in out if r['MVP_exact_pair_hits']==r['FinnGen_exact_pair_hits']==1];chosen=[r for r in good if r['both_nominal']]
print(json.dumps({'n':len(out),'valid_both':len(good),'nominal_both':len(chosen),'chosen':chosen,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()},indent=2))
