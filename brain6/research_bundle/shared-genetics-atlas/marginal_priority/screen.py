import csv,gzip,math,json,pathlib,hashlib,collections
base=pathlib.Path('/Volumes/Extreme SSD/sleep-frailty-full-autonomous-v5-20261009/frailty_paper/full_autonomous_research_v5/genetics/new_public_validation/common_support');out=pathlib.Path('outputs/shared-genetics-atlas/marginal_priority');names=['S','FI','FFS'];files=[base/(x+'.sumstats.gz') for x in names];n=1138262;m=3*n;threshold=.05/m;hits=[];seen=set();counts=collections.Counter()
with gzip.open(files[0],'rt') as a,gzip.open(files[1],'rt') as b,gzip.open(files[2],'rt') as c:
 for rows in zip(*(csv.DictReader(f,delimiter='\t') for f in (a,b,c)),strict=True):
  key=tuple(rows[0][k] for k in ('SNP','A1','A2'));assert all(tuple(t[k] for k in ('SNP','A1','A2'))==key for t in rows);assert key[0] not in seen;seen.add(key[0]);z=[float(t['Z']) for t in rows];assert all(math.isfinite(x) for x in z);pv=[math.erfc(max(abs(x)-.0005,0)/math.sqrt(2)) for x in z]; passed=[names[i] for i in range(3) if pv[i]<=threshold]
  for trait in passed:counts[trait]+=1
  if passed:hits.append(dict(SNP=key[0],A1=key[1],A2=key[2],Z_insomnia=z[0],Z_FI=z[1],Z_FFS=z[2],p_insomnia_upper=pv[0],p_FI_upper=pv[1],p_FFS_upper=pv[2],threshold_crossing='+'.join(passed)))
assert len(seen)==n
fi=pathlib.Path('/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/frailty_v1/harmonized/frailty.harmonized.tsv.gz'); wanted={x['SNP'] for x in hits};positions=collections.defaultdict(list)
with gzip.open(fi,'rt') as f:
 for x in csv.DictReader(f,delimiter='\t'):
  if x['SNP'] in wanted:positions[x['SNP']].append((x['CHR'],x['BP'],x['A1'],x['A2']))
for x in hits:
 matches=positions[x['SNP']];alleles={x['A1'],x['A2']};ok=[p for p in matches if {p[2],p[3]}==alleles]
 x['chr37'],x['bp37'],x['position_status']=(ok[0][0],ok[0][1],'UNIQUE_ALLELE_MATCH') if len(matches)==len(ok)==1 else ('','','MISSING_OR_AMBIGUOUS')
counts['either']=len(hits);counts.update(collections.Counter('pattern_'+x['threshold_crossing'] for x in hits));counts['unique_position']=sum(x['position_status']=='UNIQUE_ALLELE_MATCH' for x in hits)
cols=['SNP','chr37','bp37','position_status','A1','A2','Z_insomnia','Z_FI','Z_FFS','p_insomnia_upper','p_FI_upper','p_FFS_upper','threshold_crossing'];
with (out/'marginal_threshold_crossings.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=cols,delimiter='\t');w.writeheader();w.writerows(sorted(hits,key=lambda x:min(x['p_insomnia_upper'],x['p_FI_upper'],x['p_FFS_upper'])))
meta={'n_support':n,'n_hypotheses':m,'threshold':threshold,'counts':counts,'input_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in files},'position_source':str(fi),'position_source_bytes':fi.stat().st_size}
(out/'summary.json').write_text(json.dumps(meta,indent=2)+'\n');print(json.dumps(meta,indent=2))
