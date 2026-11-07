import csv,pathlib,json,collections
p=pathlib.Path(__file__).resolve().parent;allrows=[];summary=[];compare_cols=['chromosome','base_pair_location','effect_allele','other_allele','beta','standard_error','effect_allele_frequency','p_value','rsid','rs_id','n','Q_metric','Q_metric_pvalue','hm_coordinate_conversion','hm_code','variant_id']
for n in range(90624046,90624053):
 acc=f'GCST{n}'
 with (p/f'{acc}_extended.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
 with (p/f'{acc}_anchors.tsv').open() as f:anchors=list(csv.DictReader(f,delimiter='\t'))
 assert len(rows)==87 and len({x['selected_rsid'] for x in rows})==87
 byid={x['selected_rsid']:x for x in rows}
 for x in anchors:
  y=byid[x['rsid_requested']];assert x['status']==y['status']
  if x['status']=='UNIQUE_MATCH':
   for key in compare_cols:assert x[key]==y[key],(acc,x['rsid_requested'],key)
 for x in rows:
  x['below_published_factor_threshold']='NA' if x['status']!='UNIQUE_MATCH' else str(float(x['p_value'])<7.14e-9)
  allrows.append(x)
 hits=[x for x in rows if x['below_published_factor_threshold']=='True']
 summary.append({'accession':acc,'status_counts':dict(collections.Counter(x['status'] for x in rows)),'passing_marker_count':len(hits),'passing_marker_chromosomes':dict(collections.Counter(x['chr38'] for x in hits)),'five_anchor_comparison':'identical statuses and returned source fields'})
with (p/'all_extended.tsv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(allrows[0]),delimiter='\t');w.writeheader();w.writerows(allrows)
(p/'extended_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print('validated',len(allrows),'selected variant-factor rows and 35 separate same-source lookup comparisons')
