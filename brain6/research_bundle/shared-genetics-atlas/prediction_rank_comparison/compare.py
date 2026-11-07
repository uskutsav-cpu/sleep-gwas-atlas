import csv,collections,math,pathlib,json,hashlib
base=pathlib.Path('outputs/brain6-alphagenome-pilot/results');out=pathlib.Path('outputs/shared-genetics-atlas/prediction_rank_comparison');p=base/'jbrowse_expanded/expanded_variant_evidence.tsv';r=base/'regions.tsv'
with p.open() as f: rows=list(csv.DictReader(f,delimiter='\t'))
with r.open() as f: regions={x['candidate_locus_id']:x for x in csv.DictReader(f,delimiter='\t')}
assert len(rows)==425 and len(regions)==6 and all(x['tissue_prediction_status']=='PREDICTED' for x in rows)
groups=collections.defaultdict(list)
for x in rows:groups[x['candidate_locus_id']].append(x)
def ranks(items,fun):
 vals=[fun(x) for x in items];assert all(math.isfinite(v) for v in vals);sv=sorted(vals,reverse=True);return {x['snp']:sv.index(v)+1 for x,v in zip(items,vals)}
ranked=[]; summary=[]
for rid,items in groups.items():
 assert len(items)==int(regions[rid]['n_union_credible_variants']) and len({x['snp'] for x in items})==len(items)
 dr=ranks(items,lambda x:abs(float(x['dnase_log2_ratio_501bp'])));avi=ranks(items,lambda x:float(x['avi_raw_e5']))
 pp={}
 for trait in ('insomnia','adhd'):
  valid=regions[rid][trait+'_set_status']=='EXPLORATORY_ABF_CREDIBLE_SET'
  scored=[x for x in items if x[trait+'_conditional_abf_pip']!=''];pp[trait]=ranks(scored,lambda x:float(x[trait+'_conditional_abf_pip'])) if valid and scored else None
  dtop={x['snp'] for x in items if dr[x['snp']]<=5};atop={x['snp'] for x in items if avi[x['snp']]<=5}
  ptop={x['snp'] for x in items if pp[trait] and pp[trait].get(x['snp'],999999)<=5}
  summary.append({'region':rid,'trait':trait,'n':len(items),'n_pip_available':len(scored) if valid else '', 'trait_rank_gate':'PASSED' if valid else 'NOT_ESTIMATED','dnase_top5_n':len(dtop),'pip_top5_n':len(ptop) if valid else '','dnase_pip_top5_overlap':len(dtop&ptop) if valid else '','dnase_pip_shared':';'.join(sorted(dtop&ptop)) if valid else '','dnase_avi_top5_overlap':len(dtop&atop),'dnase_avi_shared':';'.join(sorted(dtop&atop))})
 for x in items:
  ranked.append({'region':rid,'snp':x['snp'],'grch38':x['atlas_variant'],'dnase_alt_ref_log2_501bp':x['dnase_log2_ratio_501bp'],'dnase_absolute_rank':dr[x['snp']],'avi_raw_rank':avi[x['snp']],'insomnia_pip':x['insomnia_conditional_abf_pip'],'insomnia_pip_rank':pp['insomnia'].get(x['snp'],'') if pp['insomnia'] else '', 'adhd_pip':x['adhd_conditional_abf_pip'],'adhd_pip_rank':pp['adhd'].get(x['snp'],'') if pp['adhd'] else ''})
for name,data in [('variant_ranks.tsv',ranked),('overlap_summary.tsv',summary)]:
 with (out/name).open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(data[0]),delimiter='\t');w.writeheader();w.writerows(data)
(out/'input_hashes.json').write_text(json.dumps({str(x):hashlib.sha256(x.read_bytes()).hexdigest() for x in (p,r)},indent=2)+'\n')
for x in summary:print(x)
for x in ranked:
 if x['snp'] in ('rs62264764','rs2431108','rs4916723','rs1476535'):print('candidate',x)
