import csv,hashlib,json,re
from pathlib import Path
p=Path(__file__).resolve().parent
root=Path('/Users/swethasunilkumar/Documents/Codex/2026-10-27/brain6')
tree=json.loads((p/'repository_tree.json').read_text()); commit='6f370810052d08d69ba465d205b34966fc59be2a'; assert tree['sha']==commit and tree['truncated'] is False
names=[x['path'] for x in tree['tree'] if x['type']=='blob'];assert not any(x.startswith('Data/') for x in names)
sens=(root/'work/mak_2.2_MR_frailty_sensitivity.R').read_bytes();main=(root/'work/mak_2.1_MR_frailty_selected_metabolites.R').read_bytes();a=sens.decode();b=main.decode()
assert '"nomental"' in a and 'Frailty/fi_",mod_fi_list[i],".fastGWA"' in a and '"_no_metabolomics.fastGWA"' in a and '"_random50.fastGWA"' in a
assert all(x in a for x in ['snp_col = "SNP"','beta_col = "BETA"','se_col = "SE"','eaf_col = "AF1"','effect_allele_col = "A1"','other_allele_col = "A2"','pval_col = "P"','samplesize_col = "N"'])
assert 'fi.fastGWA' in b and 'N=403,041' in b
meta=json.loads((p.parent/'figshare_v1_metadata.json').read_text()); file_names=[x['name'] for x in meta['files']]
assert 'fi.fastGWA' in file_names and not any('nomental' in x.lower() for x in file_names)
sha=lambda x:hashlib.sha256(x).hexdigest()
r={'commit':commit,'repository_tree_sha256':sha((p/'repository_tree.json').read_bytes()),'script_sha256':{'2.2_MR_frailty_sensitivity.R':sha(sens),'2.1_MR_frailty_selected_metabolites.R':sha(main)},'nomental_paths_implied_by_code':['Data/GWAS_summary_statistics_data/Frailty/fi_nomental.fastGWA','Data/GWAS_summary_statistics_data/Frailty/fi_nomental_no_metabolomics.fastGWA','Data/GWAS_summary_statistics_data/Frailty/fi_nomental_random50.fastGWA'],'public_tree_blob_count':len(names),'public_tree_data_blobs':0,'figshare_v1_filenames':len(file_names),'full_fi_present':True,'nomental_named_files':False,'scope':'Tree and one versioned dataset inventory checked; does not establish absence in all repositories or rights to unpublished statistics.'}
(p/'verification.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
