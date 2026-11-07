import csv,json,hashlib
from pathlib import Path
root=Path(__file__).resolve().parents[3]; out=Path(__file__).parent
paths={'lee':root/'shared-genetics-atlas/variant_context/crossdisorder_mpra/published_table/candidate_comparison.tsv','pred':root/'brain6-alphagenome-pilot/results/jbrowse_expanded/expanded_variant_evidence.tsv','window':root/'shared-genetics-atlas/alphagenome_window_sensitivity/window_scores.tsv'}
def read(p):
 with p.open(newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
a=read(paths['lee']); b=read(paths['pred']); c=read(paths['window']); by={x['snp']:x for x in b}; ww={x['snp']:x for x in c}; rows=[]
for x in a:
 if x['allelic_status']!='unique_match' or x['element_status']!='unique_match': continue
 v=by[x['rsid']]; w=ww[x['rsid']]; assert x['grch38']==v['atlas_variant']==w['variant']; assert v['prediction_window_sha256']==w['input_sha256']
 rows.append({'rsid':x['rsid'],'region':v['candidate_locus_id'],'variant38':v['atlas_variant'],'lee_allelic_class':x['allelic_class'],'lee_allelic_fdr':x['allelic_fdr'],'lee_element_class':x['element_class'],'lee_element_fdr':x['element_fdr'],'dnase_101':w['log2ratio_101_pc1'],'dnase_501':w['log2ratio_501_pc1'],'dnase_2001':w['log2ratio_2001_pc1'],'adult_rna_501':v['adult_rna_log2_ratio_501bp'],'embryonic_plus_501':v['embryonic_plus_rna_log2_ratio_501bp'],'embryonic_minus_501':v['embryonic_minus_rna_log2_ratio_501bp'],'insomnia_pip':v['insomnia_conditional_abf_pip'],'adhd_pip':v['adhd_conditional_abf_pip']})
assert len(rows)==27 and len({x['rsid'] for x in rows})==27
with (out/'lee_27_crosswalk.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=rows[0],delimiter='\t');w.writeheader();w.writerows(rows)
(out/'lee_27_receipt.json').write_text(json.dumps({'rows':len(rows),'inputs':{k:hashlib.sha256(v.read_bytes()).hexdigest() for k,v in paths.items()}},indent=2)+'\n')
for x in rows:
 if x['lee_allelic_class']=='emVar': print(x)
