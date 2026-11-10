import csv,json,glob,datetime,hashlib
from pathlib import Path
out=Path('outputs'); queue=list(csv.DictReader(open(out/'stage5_score_queue.csv')))
found={}
for file in glob.glob('work/stage5_indexed_results/*.json'):
 x=json.load(open(file)); found.update(x['found'])
snv=[r for r in queue if len(r['reference_allele'])==len(r['alternate_allele'])==1]
missing=[r['variant_key'] for r in snv if r['variant_key'] not in found]
assert len(found)==len(snv)==460 and not missing,(len(found),len(snv),missing)
url='https://grr.iossifovlab.com/hg38/scores/AVI/alphagenome_variant_impact_score_snvs.tsv.gz'
with open(out/'stage5_all_snv_atlas_avi_precise.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=['variant_GRCh38_ref_alt','input_classes','input_labels','avi_raw','avi_phred','source_url','source_file_md5','index_md5','note']);w.writeheader()
 for r in snv:
  v=found[r['variant_key']];w.writerow(dict(variant_GRCh38_ref_alt=r['variant_key'],input_classes=r['input_classes'],input_labels=r['input_labels'],avi_raw=v['raw_score'],avi_phred=v['avi_phred'],source_url=url,source_file_md5='d7753d6a3aba87a0b90ecb9df3b789e3',index_md5='911ee04b6e8fba1882465221b5a14f2c',note='General predicted variant impact; not disease risk, not tissue-specific; published decimal precision'))
web={}
for r in csv.DictReader(open(out/'stage5_23_markers_atlas_avi.csv')): web[r['variant_GRCh38_ref_alt']]=float(r['avi_phred_displayed_rounded'])
for r in csv.DictReader(open(out/'stage5_priority_credible_snv_atlas_avi.csv')): web[r['variant_GRCh38_ref_alt']]=float(r['avi_phred_displayed_rounded'])
for r in csv.DictReader(open(out/'stage5_avi_matched_control_comparison.csv')): web[r['control_variant']]=float(r['control_avi_phred_displayed'])
for line in open('work/stage5_atlas_additional_avi.jsonl'):
 r=json.loads(line)
 if r['status']=='ok': web[r['variant']]=float(r['avi_phred_displayed'])
comparisons=[{'variant':k,'web_display':v,'indexed_phred':float(found[k]['avi_phred']),'abs_difference':round(abs(v-float(found[k]['avi_phred'])),5)} for k,v in web.items() if k in found]
with open(out/'stage5_avi_web_index_crosscheck.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=['variant','web_display','indexed_phred','abs_difference']);w.writeheader();w.writerows(comparisons)
old=list(csv.DictReader(open(out/'stage5_preliminary_functional_atlas.csv')))
fields=list(old[0])+['atlas_avi_raw_indexed','atlas_avi_phred_indexed','indexed_avi_source']
with open(out/'stage5_expanded_functional_atlas.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
 for r in old:
  v=found.get(r['variant_GRCh38_ref_alt']);r['atlas_avi_raw_indexed']=v['raw_score'] if v else '';r['atlas_avi_phred_indexed']=v['avi_phred'] if v else '';r['indexed_avi_source']=url if v else '';w.writerow(r)
audit={'retrieval_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'total_queue':len(queue),'snvs':len(snv),'snvs_found':len(found),'indels_not_in_source':len(queue)-len(snv),'groups':len(glob.glob('work/stage5_indexed_results/*.json')),'web_compared':len(comparisons),'web_discrepancies_over_0_05001':[r for r in comparisons if r['abs_difference']>0.05001],'index_md5_verified':'911ee04b6e8fba1882465221b5a14f2c','file_md5_reported_by_source_not_independently_verified':'d7753d6a3aba87a0b90ecb9df3b789e3','note':'Indexed ranges only; score is generalized variant impact and not risk, cell type, mechanism or causal evidence.'}
json.dump(audit,open(out/'stage5_indexed_avi_audit.json','w'),indent=2)
print(json.dumps(audit,indent=2))
