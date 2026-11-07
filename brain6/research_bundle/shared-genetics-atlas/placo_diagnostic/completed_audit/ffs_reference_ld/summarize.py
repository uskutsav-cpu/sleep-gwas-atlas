import csv,json,pathlib
p=pathlib.Path(__file__).resolve().parent
summary=[]
for chrom in ('5','11','16','17','22'):
 d=json.loads((p/f'chr{chrom}_receipt.json').read_text()); pairs=list(csv.DictReader((p/f'chr{chrom}_pairs.tsv').open(),delimiter='\t')); anchor=d['anchor']
 ffsids={r['variant_1'] for r in pairs if r['pair_1']=='FFS'} | {r['variant_2'] for r in pairs if r['pair_2']=='FFS'}
 valid={x for x in ffsids if d['variants'][x]['matches']==1}; components={x:{x} for x in valid}
 for r in pairs:
  if r['status']=='OK' and r['pair_1']==r['pair_2']=='FFS' and float(r['r2'])>=.8:
   union=components[r['variant_1']]|components[r['variant_2']]
   for x in union:components[x]=union
 components=sorted({tuple(sorted(v)) for v in components.values()},key=lambda x:(-len(x),x))
 a=[float(r['r2']) for r in pairs if r['status']=='OK' and r['pair_1']==r['pair_2']=='FFS' and anchor in (r['variant_1'],r['variant_2'])]
 summary.append({'chromosome':chrom,'selected_ffs':len(ffsids),'valid_reference_records':len(valid),'anchor':anchor,'anchor_other_n':len(a),'anchor_r2_min':min(a) if a else None,'anchor_r2_max':max(a) if a else None,'r2_08_graph_components':[list(x) for x in components],'missing_reference':sorted(ffsids-valid),'pair_counts':d['pair_counts']})
(p/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
