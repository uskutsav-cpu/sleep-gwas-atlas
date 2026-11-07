import csv,hashlib,json,math,pathlib,pysam
p=pathlib.Path(__file__).resolve().parent;root=p.parents[1];ctx=root/'variant_context/selected_variant_gene_context.tsv';url='https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90568001-GCST90569000/GCST90568440/harmonised/GCST90568440.h.tsv.gz'
selected=['rs62264764','rs77960','rs1476535'];expected={'rs62264764':('3',117920728,'G','A'),'rs77960':('5',104628884,'G','A'),'rs1476535':('7',114430980,'C','T')};fields=['chromosome','base_pair_location','effect_allele','other_allele','beta','standard_error','effect_allele_frequency','p_value','rsid','hm_coordinate_conversion','hm_code','variant_id']
with ctx.open(newline='') as h:context={r['snp']:r for r in csv.DictReader(h,delimiter='\t') if r['snp'] in selected}
assert len(context)==3;t=pysam.TabixFile(url);out=[];receipts=[]
for snp in selected:
 ch,pos,ref,alt=expected[snp];c=context[snp];assert c['grch38_ref_alt']==f'chr{ch}:{pos}:{ref}>{alt}' and {ref,alt} not in ({'A','T'},{'C','G'})
 raw=list(t.fetch(ch,pos-1,pos));parsed=[]
 for line in raw:
  v=line.split('\t');assert len(v)==len(fields);parsed.append(dict(zip(fields,v)))
 matches=[z for z in parsed if z['chromosome']==ch and z['base_pair_location']==str(pos) and z['rsid']==snp and {z['effect_allele'],z['other_allele']}=={ref,alt}]
 row={'variant':snp,'coordinate_grch38':f'{ch}:{pos}:{ref}>{alt}','effect_allele':c['source_A1'],'discovery_adhd_beta':c['adhd_beta_A1'],'status':'NOT_ESTIMABLE','symptom_beta_aligned':'','symptom_se':'','symptom_p':'','symptom_p_times3':'','direction_agreement':''}
 if len(matches)==1:
  q=matches[0]
  try:
   beta,se,pv,f=map(float,(q['beta'],q['standard_error'],q['p_value'],q['effect_allele_frequency']))
   assert all(map(math.isfinite,(beta,se,pv,f))) and se>0 and 0<=pv<=1 and 0<=f<=1
   aligned=beta*(1 if q['effect_allele']==c['source_A1'] else -1);orig=float(c['adhd_beta_A1'])
   row.update(status='MATCHED',symptom_beta_aligned=aligned,symptom_se=se,symptom_p=pv,symptom_p_times3=min(1,3*pv),direction_agreement='YES' if aligned*orig>0 else ('NO' if aligned*orig<0 else 'ZERO'))
  except (ValueError,AssertionError): row['status']='INVALID_STATISTICS'
 elif len(matches)>1:row['status']='AMBIGUOUS'
 out.append(row);receipts.append({'variant':snp,'region':f'{ch}:{pos}-{pos}','raw_rows':raw,'matches':len(matches),'raw_rows_sha256':hashlib.sha256(('\n'.join(raw)+'\n').encode()).hexdigest()})
with (p/'selected_variants.tsv').open('w',newline='') as h:w=csv.DictWriter(h,fieldnames=out[0].keys(),delimiter='\t');w.writeheader();w.writerows(out)
receipt={'url':url,'metadata_sha256':hashlib.sha256((p/'metadata.yaml').read_bytes()).hexdigest(),'context_sha256':hashlib.sha256(ctx.read_bytes()).hexdigest(),'query_records':receipts,'note':'Indexed remote queries; full compressed GWAS and tabix index not hashed. Individual-cohort overlap not established.'};(p/'query_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(out,indent=2))
