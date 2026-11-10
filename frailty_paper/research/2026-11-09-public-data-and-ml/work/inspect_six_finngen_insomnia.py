import csv,gzip,hashlib,pathlib,collections,json,math,re
p=pathlib.Path('/Volumes/Extreme SSD/sleep-frailty-full-autonomous-v5-20261009/frailty_paper/full_autonomous_research_v5/genetics/new_public_validation/raw/finngen_R13_F5_INSOMNIA.gz'); sh=hashlib.sha256()
with p.open('rb') as f:
 for b in iter(lambda:f.read(8*1024*1024),b''):sh.update(b)
assert sh.hexdigest()=='353129dc8252461a2a97087ce6de39548d0b0bf1e488a6ef36da75e9f4690049'
r=list(csv.DictReader(open('outputs/song124_mvp_hfrs_joint_screen.csv')));hits={x['variant']:x for x in r if x.get('joint_bonferroni_124') and float(x['joint_bonferroni_124'])<.05};assert len(hits)==6
found=collections.defaultdict(list);n=0
with gzip.open(p,'rt') as f:
 h=f.readline().rstrip('\r\n').split('\t');idx=h.index('rsids')
 for line in f:
  n+=1;v=line.rstrip('\r\n').split('\t'); ids=set(re.split(r'[,;\s]+',v[idx])) & hits.keys()
  for rs in ids:found[rs].append(dict(zip(h,v)))
out=[]
for rs,c in hits.items():
 same=[x for x in found[rs] if x['#chrom']==c['GRCh38_chr'] and x['pos']==c['GRCh38_pos'] and {x['ref'],x['alt']}=={c['HFRS_effect_allele'],c['HFRS_other_allele']}]
 r={'variant':rs,'GRCh38_chr':c['GRCh38_chr'],'GRCh38_pos':c['GRCh38_pos'],'hfrs_effect_allele':c['HFRS_effect_allele'],'hfrs_other_allele':c['HFRS_other_allele'],'id_rows':len(found[rs]),'exact_position_allele_rows':len(same)}
 if len(same)==1:
  x=same[0]; b=float(x['beta'])*(1 if x['alt']==c['HFRS_effect_allele'] else -1); se=float(x['sebeta']);r.update(r13_beta_aligned=b,r13_se=se,r13_ci_low=b-1.95996398454*se,r13_ci_high=b+1.95996398454*se,r13_p=x['pval'],r13_alt=x['alt'],r13_ref=x['ref'],r13_alt_frequency=x['af_alt'],source_nearest_genes=x['nearest_genes'])
 out.append(r)
keys=list(dict.fromkeys(k for r in out for k in r)); dest=pathlib.Path('outputs/six_joint_hits_finngen_insomnia_followup.csv')
with dest.open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(out)
print('source_rows',n,'sha256',sh.hexdigest(),'output_sha256',hashlib.sha256(dest.read_bytes()).hexdigest());print(json.dumps(out,indent=2))
