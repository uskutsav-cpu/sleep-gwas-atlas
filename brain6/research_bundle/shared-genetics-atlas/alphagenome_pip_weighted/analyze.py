import csv,collections,pathlib,math,hashlib,json
p=pathlib.Path(__file__).resolve().parent;atlas=p.parent;variant=atlas.parent/'brain6-alphagenome-pilot/results/jbrowse_expanded/expanded_variant_evidence.tsv';regions=atlas.parent/'brain6-alphagenome-pilot/results/regions.tsv'
with variant.open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
with regions.open() as f:metadata={x['candidate_locus_id']:x for x in csv.DictReader(f,delimiter='\t')}
assert len(rows)==425 and len(metadata)==6
by=collections.defaultdict(list)
for x in rows:by[x['candidate_locus_id']].append(x)
summary=[];contributions=[]
for region,xs in sorted(by.items()):
 assert len(xs)==int(metadata[region]['n_union_credible_variants'])
 for trait in ('insomnia','adhd'):
  if metadata[region][f'{trait}_set_status']!='EXPLORATORY_ABF_CREDIBLE_SET':continue
  vals=[];missing=0
  for x in xs:
   raw=x[f'{trait}_conditional_abf_pip']
   if not raw.strip():
    missing+=1;continue
   pp=float(raw);effect=float(x['dnase_log2_ratio_501bp']);assert math.isfinite(pp) and 0<=pp<=1 and math.isfinite(effect);vals.append((x['snp'],pp,effect))
  mass=math.fsum(x[1] for x in vals);numer=math.fsum(pp*abs(effect) for _,pp,effect in vals);assert 0<mass<=1.000001 and numer>0
  assert abs(math.fsum(pp*abs(effect)/numer for _,pp,effect in vals)-1)<1e-12
  ranked=sorted(vals,key=lambda x:(-x[1]*abs(x[2]),x[0]));winner,wp,we=ranked[0]
  summary.append({'region':region,'trait':trait,'selected_union_variants':len(xs),'pip_available_variants':len(vals),'pip_missing_variants':missing,'selected_pip_mass':mass,'unweighted_mean_abs':sum(abs(x[2]) for x in vals)/len(vals),'pip_weighted_mean_abs':numer/mass,'largest_contribution_variant':winner,'largest_contribution_pip':wp,'largest_contribution_abs':abs(we),'largest_contribution_share':wp*abs(we)/numer})
  for snp,pp,effect in vals:contributions.append({'region':region,'trait':trait,'snp':snp,'pip':pp,'dnase_signed':effect,'abs_dnase':abs(effect),'normalized_pip':pp/mass,'weighted_contribution':pp*abs(effect)/mass,'share_of_weighted_sum':pp*abs(effect)/numer})
for name,rows in [('region_summary.tsv',summary),('variant_contributions.tsv',contributions)]:
 with (p/name).open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
(p/'receipt.json').write_text(json.dumps({'input_sha256':{str(x):sha(x) for x in (variant,regions)},'variants':425,'region_trait_pairs':len(summary),'method':'PIP-weighted absolute predicted local adult DNase ratio, conditional normalization over selected union only','restriction':'No AlphaGenome outputs used for training; exploratory scientific interpretation only.'},indent=2)+'\n');print('analyzed',len(summary),'admitted region–trait pairs')
