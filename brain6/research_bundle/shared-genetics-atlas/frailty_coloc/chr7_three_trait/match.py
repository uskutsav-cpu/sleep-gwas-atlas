import csv,collections,hashlib,json,math,pathlib
p=pathlib.Path(__file__).resolve().parent;root=p.parents[1];src=pathlib.Path('/Volumes/Extreme SSD/brain6-work/exploratory-finemap-coloc-v2/inputs');out=pathlib.Path('/Volumes/Brain6Work/chr7-three-trait-input.tsv')
def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
def read(f):
 with f.open(newline='') as h:return list(csv.DictReader(h,delimiter='\t'))
fi_path=root/'frailty_coloc/matched.tsv'; fi=[r for r in read(fi_path) if r['chr']=='7'];prov=json.loads((root.parent/'brain6-alphagenome-pilot/results/provenance.json').read_text())['input_sha256'];indices={};sources={};counts={}
for trait in ('insomnia','adhd'):
 name='placo_insomnia__adhd_chr7_113571035_114571035__'+trait+'.tsv';f=src/name;assert sha(f)==prov[name];rows=read(f);idx=collections.defaultdict(list)
 for r in rows:idx[r['SNP']].append(r)
 indices[trait]=idx;sources[trait]={'path':str(f),'sha256':sha(f),'n':len(rows)}
assert len({r['rsid'] for r in fi})==len(fi)
res=[];excluded=collections.Counter()
for a in fi:
 matches={};reason=None
 for trait in ('insomnia','adhd'):
  rs=indices[trait].get(a['rsid'],[])
  if len(rs)!=1:reason=trait+'_not_unique';break
  b=rs[0]
  if b['CHR']!='7' or b['POS']!=a['bp37']:reason=trait+'_coordinate';break
  alleles={b['A1'],b['A2']}
  if len(alleles)!=2 or not alleles.issubset(set('ACGT')) or alleles in ({'A','T'},{'C','G'}) or alleles!={a['effect_allele'],a['other_allele']}:reason=trait+'_alleles';break
  try:
   beta,se=float(b['BETA']),float(b['SE'])
   if not math.isfinite(beta) or not math.isfinite(se) or se<=0:raise ValueError()
  except ValueError:reason=trait+'_statistics';break
  matches[trait]=(beta*(1 if b['A1']==a['effect_allele'] else -1),se)
 if reason:excluded[reason]+=1;continue
 beta,se=float(a['fi_beta']),float(a['fi_se'])
 if not math.isfinite(beta) or not math.isfinite(se) or se<=0:excluded['fi_statistics']+=1;continue
 res.append({'rsid':a['rsid'],'pos37':a['bp37'],'effect_allele':a['effect_allele'],'other_allele':a['other_allele'],'fi_beta':beta,'fi_se':se,'insomnia_beta':matches['insomnia'][0],'insomnia_se':matches['insomnia'][1],'adhd_beta':matches['adhd'][0],'adhd_se':matches['adhd'][1]})
assert len(res)>=500
with out.open('w',newline='') as h:w=csv.DictWriter(h,fieldnames=res[0].keys(),delimiter='\t');w.writeheader();w.writerows(res)
receipt={'fi_input_path':str(fi_path),'fi_input_sha256':sha(fi_path),'fi_chr7_rows':len(fi),'brain6_sources':sources,'excluded_sequential':dict(excluded),'matched':len(res),'external_analysis_input':str(out),'external_analysis_input_sha256':sha(out),'note':'Exact rsID/position/allele match and orientation only; not an independent source-frequency, model or overlap validation.'}
(p/'matching_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
