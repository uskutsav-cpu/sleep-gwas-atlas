import csv,json,urllib.request,urllib.parse,time,hashlib,pathlib,collections,subprocess
base=pathlib.Path('/Users/swethasunilkumar/Documents/Codex/2026-10-27/brain6'); src=base/'outputs/brain6-alphagenome-pilot/results/jbrowse_expanded/expanded_variant_evidence.tsv';out=base/'outputs/shared-genetics-atlas/variant_context'; raw=out/'ensembl';raw.mkdir(parents=True,exist_ok=True)
with src.open() as f: rows=list(csv.DictReader(f,delimiter='\t'))
groups=collections.defaultdict(list)
for r in rows:groups[r['candidate_locus_id']].append(r)
selected={}
for reg,items in groups.items():
 for label,col,absolute in [('max_abs_adult_dnase','dnase_log2_ratio_501bp',True),('max_avi','avi_phred_lower',False),('max_insomnia_pip','insomnia_conditional_abf_pip',False),('max_adhd_pip','adhd_conditional_abf_pip',False)]:
  best=max(items,key=lambda a:abs(float(a[col] or 0)) if absolute else float(a[col] or 0)); key=(reg,best['snp']); selected.setdefault(key,{'row':best,'reasons':[]})['reasons'].append(label)
results=[]
for (reg,snp),v in selected.items():
 r=v['row']; chrom=r['chromosome'].removeprefix('chr'); pos=int(r['pos_grch38_1based']);url=f'https://rest.ensembl.org/overlap/region/human/{chrom}:{pos-100000}-{pos+100000}?feature=gene;content-type=application/json'; path=raw/(snp+'.json')
 try:
  if not path.exists():
   data=subprocess.run(['curl','--fail','--silent','--show-error','--max-time','35','-H','Accept: application/json',url],check=True,capture_output=True).stdout;parsed=json.loads(data);path.write_bytes(data);time.sleep(.2)
  else:data=path.read_bytes();parsed=json.loads(data)
  if not isinstance(parsed,list):raise ValueError(str(parsed)[:100])
  genes=[]
  for g in parsed:
   if g.get('start') is None or g.get('end') is None:continue
   distance=max(g['start']-pos,pos-g['end'],0)
   genes.append({'name':g.get('external_name',g.get('id','')),'ensembl_id':g.get('id'),'biotype':g.get('biotype'),'distance_to_gene_body_bp':distance,'start':g['start'],'end':g['end']})
  genes.sort(key=lambda g:(g['distance_to_gene_body_bp'],g['name']))
  coding=[g for g in genes if g['biotype']=='protein_coding'];nearest=genes[0] if genes else {};nearest_coding=coding[0] if coding else {}
  summary={'status':'ok','gene_count':len(genes),'overlapping_genes':';'.join(g['name']+'|'+g['ensembl_id'] for g in genes if g['distance_to_gene_body_bp']==0),'nearest_gene':nearest.get('name',''),'nearest_gene_distance_bp':nearest.get('distance_to_gene_body_bp',''),'nearest_protein_coding_gene':';'.join(g['name']+'|'+g['ensembl_id'] for g in coding if g['distance_to_gene_body_bp']==nearest_coding.get('distance_to_gene_body_bp')),'nearest_coding_distance_bp':nearest_coding.get('distance_to_gene_body_bp',''),'raw_response_sha256':hashlib.sha256(data).hexdigest()}
 except Exception as e:summary={'status':'ERROR '+str(e)}
 results.append({'region':reg,'snp':snp,'selection_reasons':';'.join(v['reasons']),'grch38_ref_alt':r['atlas_variant'],'source_A1':r['source_a1'],'insomnia_beta_A1':r['insomnia_beta_a1'],'adhd_beta_A1':r['adhd_beta_a1'],'insomnia_pip':r['insomnia_conditional_abf_pip'],'adhd_pip':r['adhd_conditional_abf_pip'],'insomnia_credible_set_member':r['insomnia_credible_set_member'],'adhd_credible_set_member':r['adhd_credible_set_member'],'avi_phred':r['avi_phred_lower'],'adult_dnase_log2ratio':r['dnase_log2_ratio_501bp'],'query_url':url,**summary})
 print(snp,summary.get('status'),summary.get('overlapping_genes',''),flush=True)
keys=[]
for r in results:
 for k in r:
  if k not in keys:keys.append(k)
with (out/'selected_variant_gene_context.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,keys,delimiter='\t');w.writeheader();w.writerows(results)
print('Selected',len(results),'errors',sum(not r['status']=='ok' for r in results))
