import csv,gzip,hashlib,pathlib,collections,json,math,re
p=pathlib.Path('/Volumes/Extreme SSD/sleep-frailty-full-autonomous-v5-20261009/frailty_paper/full_autonomous_research_v5/genetics/new_public_validation/raw/finngen_R13_F5_INSOMNIA.gz')
sha=hashlib.sha256()
with p.open('rb') as f:
 for b in iter(lambda:f.read(8*1024*1024),b''):sha.update(b)
assert sha.hexdigest()=='353129dc8252461a2a97087ce6de39548d0b0bf1e488a6ef36da75e9f4690049'
rows=list(csv.DictReader(open('outputs/expanded_22_candidate_risk_gene_table.csv')));w={x['variant']:x for x in rows};coords={(x['GRCh38_chr'],x['GRCh38_pos']):x['variant'] for x in rows};assert len(coords)==22
found=collections.defaultdict(list);n=0
with gzip.open(p,'rt') as f:
 h=f.readline().rstrip('\r\n').split('\t');i={k:h.index(k) for k in ('#chrom','pos','rsids','ref','alt')}
 for line in f:
  n+=1;v=line.rstrip('\r\n').split('\t');rs=coords.get((v[i['#chrom']],v[i['pos']]))
  if rs and rs in re.split(r'[,;\s]+',v[i['rsids']]):found[rs].append(dict(zip(h,v)))
print('Scanned',n,'records',flush=True)
out=[]
for x in rows:
 rs=x['variant'];same=[y for y in found[rs] if {y['ref'],y['alt']}=={x['counted_allele'],x['other_allele']}];d={'variant':rs,'GRCh38_chr':x['GRCh38_chr'],'GRCh38_pos':x['GRCh38_pos'],'counted_allele':x['counted_allele'],'other_allele':x['other_allele'],'nearest_gene_body':x['nearest_gene_body_positional'],'id_position_hits':len(found[rs]),'exact_allele_hits':len(same)}
 if len(same)==1:
  z=same[0];beta=float(z['beta'])*(1 if z['alt']==x['counted_allele'] else -1);se=float(z['sebeta']);mvp=float(x['MVP_insomnia_log_OR_aligned_from_rounded_OR']);d.update(r13_beta=beta,r13_se=se,r13_ci_low=beta-1.95996398454*se,r13_ci_high=beta+1.95996398454*se,r13_p=float(z['pval']),r13_ref=z['ref'],r13_alt=z['alt'],r13_alt_frequency=z['af_alt'],same_sign_as_mvp=beta*mvp>0,same_sign_all_six=all(float(x[k])*mvp>0 for k in ('FI_beta_aligned','FFS_beta_aligned','UKB_insomnia_beta_aligned','HFRS_beta_aligned')) and beta*mvp>0,all_six_nominal=all(float(x[k])<.05 for k in ('FI_p','FFS_p','UKB_insomnia_p','MVP_insomnia_p','HFRS_p')) and float(z['pval'])<.05)
 out.append(d)
dest=pathlib.Path('outputs/22_finngen_insomnia_followup.csv');keys=list(dict.fromkeys(k for d in out for k in d))
with dest.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(out)
matched=[x for x in out if x['exact_allele_hits']==1];res={'source_sha256':sha.hexdigest(),'source_rows_scanned':n,'unique_exact_matches':len(matched),'nominal_r13':sum(x['r13_p']<.05 for x in matched),'r13_same_sign_as_mvp':sum(x['same_sign_as_mvp'] for x in matched),'same_sign_and_nominal_all_six':[x['variant'] for x in matched if x['all_six_nominal'] and x['same_sign_all_six']],'output_sha256':hashlib.sha256(dest.read_bytes()).hexdigest()};pathlib.Path('outputs/22_finngen_insomnia_followup.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res,indent=2))
