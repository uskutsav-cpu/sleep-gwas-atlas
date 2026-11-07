import pathlib,hashlib,collections,csv,zlib,io,json,math,statistics,pysam
p=pathlib.Path(__file__).resolve().parent;prev=p.parent/'zenodo_prefix_audit'
panel=pathlib.Path('/Users/swethasunilkumar/Documents/Codex/2026-10-27/loo/work/1kg_phase3_panel.tsv')
with panel.open() as f:eur={r['sample'] for r in csv.DictReader(f,delimiter='\t') if r['super_pop']=='EUR'}
assert len(eur)==503
with (prev/'selected_rows.tsv').open() as f:selected=list(csv.DictReader(f,delimiter='\t'))
assert len(selected)==20
configs={'Exhaustion':('exhaustion_first_131072_bytes.bin',47139,349680),'Low_Grip_strength':('Low_Grip_strength_first_131072_bytes.bin',56262,350718),'Low_Physical_activity':('Low_Physical_activity_first_131072_bytes.bin',34406,371518),'Slow_Walking_speed':('Slow_Walking_speed_first_131072_bytes.bin',30846,375531),'Weight_loss':('Weight_loss_first_131072_bytes.bin',61378,340745)}
rows_by={};summ=[];sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest();hashes={}
for name,(filename,cases,controls) in configs.items():
 raw=(prev/filename).read_bytes();hashes[name]=hashlib.sha256(raw).hexdigest();s=zlib.decompressobj(31).decompress(raw).decode();rows=list(csv.DictReader(io.StringIO(s[:s.rfind('\n')]),delimiter='\t'));idx=collections.defaultdict(list);ratio=[];invalid=0
 for x in rows:
  idx[x['ID']].append(x)
  try:
   f=float(x['A1FREQ']);se=float(x['SE']);assert math.isfinite(f) and math.isfinite(se) and 0<f<1 and se>0
   ratio.append(se*se*2*f*(1-f)*cases*controls/(cases+controls))
  except (ValueError,AssertionError):invalid+=1
 rows_by[name]=idx;values=sorted(ratio);assert values
 def quantile(frac):return values[round(frac*(len(values)-1))]
 summ.append({'component':name,'prefix_rows':len(rows),'eligible_rows':len(values),'excluded_rows':invalid,'cases':cases,'controls':controls,'ratio_p05':quantile(.05),'ratio_median':quantile(.5),'ratio_p95':quantile(.95)})
url='https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz';vcf=pysam.VariantFile(url);out=[]
for x in selected:
 pos=int(x['GENPOS']);vs=[v for v in vcf.fetch('1',pos-1,pos) if v.pos==pos and v.id and x['ID'] in v.id.split(';') and {x['ALLELE0'],x['ALLELE1']}.issubset(v.alleles)];assert len(vs)==1
 v=vs[0]
 for name in configs:
  candidates=[r for r in rows_by[name][x['ID']] if r['CHROM']==x['CHROM'] and r['GENPOS']==x['GENPOS'] and {r['ALLELE0'],r['ALLELE1']}=={x['ALLELE0'],x['ALLELE1']}]
  result={'component':name,'rsid':x['ID'],'pos_grch37':pos,'source_matches':len(candidates)}
  if len(candidates)==1:
   y=candidates[0];allele=y['ALLELE1'];a=v.alleles.index(allele);b=v.alleles.index(y['ALLELE0']);count=0;den=0;eligible=0
   for sample in eur & set(vcf.header.samples):
    g=v.samples[sample].get('GT')
    if g and len(g)==2 and all(k in (a,b) for k in g):count+=sum(k==a for k in g);den+=2;eligible+=1
   assert den>0
   f=float(y['A1FREQ']);fref=count/den;result.update(status='MATCH',counted_allele=allele,people=eligible,source_frequency=f,reference_frequency=fref,difference=f-fref)
  else:result['status']='ABSENT_OR_NONUNIQUE'
  out.append(result)
fields=['component','rsid','pos_grch37','source_matches','status','counted_allele','people','source_frequency','reference_frequency','difference']
for filename,data,keys in [('frequency_comparisons.tsv',out,fields),('se_scale_summary.tsv',summ,list(summ[0]))]:
 with (p/filename).open('w',newline='') as f:w=csv.DictWriter(f,keys,delimiter='\t');w.writeheader();w.writerows(data)
receipt={'panel_sha256':sha(panel),'prefix_hashes':hashes,'source_selected_sha256':sha(prev/'selected_rows.tsv'),'reference_vcf':url,'frequency_output_sha256':sha(p/'frequency_comparisons.tsv'),'se_output_sha256':sha(p/'se_scale_summary.tsv'),'frequency_counts':dict(collections.Counter(x['status'] for x in out)),'caveat':'post-selected prefixes and reference subsets, not study-matched model or phenotype validation'}
(p/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
