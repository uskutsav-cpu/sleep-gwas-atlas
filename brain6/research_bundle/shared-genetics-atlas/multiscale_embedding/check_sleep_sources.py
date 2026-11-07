import csv,hashlib,json,pathlib
root=pathlib.Path('outputs/shared-genetics-atlas');out=root/'multiscale_embedding';repo=pathlib.Path('work/sleep-gwas-atlas');paths={'brain6':root/'trait_embedding/source_brain6_72.tsv','core_traits':repo/'results/atlas/traits.tsv','core_fi':repo/'frailty_paper/analysis/frozen_atlas_frailty_global_rg.tsv','fi_ffs':root/'evidence/ALL12_FI_FFS_CONTRAST_RECOMPUTED.tsv','frozen_manifest':repo/'frailty_paper/manifests/frozen_atlas_frailty_manifest.json'}
def read(p):
 with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
m=json.loads(paths['frozen_manifest'].read_text());assert sha(paths['core_traits'])==m['source_output_sha256']['results/atlas/traits.tsv'];assert sha(paths['core_fi'])==m['extracted_table_sha256']
b=read(paths['brain6']);core={x['trait_id']:x for x in read(paths['core_traits'])};fi={x['sleep_trait']:x for x in read(paths['core_fi'])};paired={x['sleep_trait']:x for x in read(paths['fi_ffs'])};seen={}
for x in b:
 t=x['sleep_trait'];v=(x['source_GWAS_ID_sleep'],x['sleep_sample_size']);assert t not in seen or seen[t]==v;seen[t]=v
assert len(seen)==12 and set(seen)==set(fi)==set(paired);rows=[]
for t,(sid,n) in seen.items():
 c=core[t];a=fi[t];p=paired[t];rows.append({'sleep_trait':t,'brain6_source_id':sid,'core_source_id':c['source_id'],'source_id_match':sid==c['source_id'],'core_dataset_version':c['dataset_version'],'brain6_sleep_n':n,'core_sleep_n':c['n_total'],'n_match':n==c['n_total'],'frozen_core_fi_rg':a['global_rg'],'matched_support_fi_rg':p['FI_rg'],'paired_support_snps':p['support_SNPs'],'conditional_minus_original_fi_rg':float(p['FI_rg'])-float(a['global_rg']),'paired_sleep_input_bytes_identity':'NOT_ESTABLISHED'})
with (out/'sleep_source_crosswalk.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
summary={'n_traits':len(rows),'source_id_matches':sum(x['source_id_match'] for x in rows),'sample_size_matches':sum(x['n_match'] for x in rows),'core_traits_sha256_verified_against_frozen_frailty_manifest':True,'core_fi_sha256_verified_against_frozen_frailty_manifest':True,'source_file_sha256':{k:sha(v) for k,v in paths.items()},'input_byte_identity_for_archived_paired_fit':'NOT_ESTABLISHED'}
(out/'sleep_source_crosswalk_receipt.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
