import csv,hashlib,json
from pathlib import Path
root=Path(__file__).resolve().parents[3]
inputs={'assays':root/'shared-genetics-atlas/variant_context/crossdisorder_mpra/published_table/cross_assay_published.tsv','predictions':root/'brain6-alphagenome-pilot/results/jbrowse_expanded/expanded_variant_evidence.tsv','windows':root/'shared-genetics-atlas/alphagenome_window_sensitivity/window_scores.tsv'}
def read(p):
 with p.open(newline='') as f: return list(csv.DictReader(f,delimiter='\t'))
a,p,w=(read(inputs[k]) for k in ('assays','predictions','windows'))
assert len(a)==4 and len(p)==425 and len(w)==425
pby={x['snp']:x for x in p}; wby={x['snp']:x for x in w}; assert len(pby)==len(wby)==425
out=[]
for x in a:
 v=pby[x['rsid']]; t=wby[x['rsid']]; assert v['atlas_variant']==t['variant'] and v['prediction_window_sha256']==t['input_sha256']
 region=v['candidate_locus_id']; local=[z for z in p if z['candidate_locus_id']==region]
 val=float(v['dnase_log2_ratio_501bp']); rank=1+sum(abs(float(z['dnase_log2_ratio_501bp']))>abs(val) for z in local)
 out.append({'rsid':x['rsid'],'region':region,'variant37':x['variant37'],'variant38':v['atlas_variant'],'region_n':len(local),'dnase_501_rank_absolute':rank,'dnase_101':t['log2ratio_101_pc1'],'dnase_501':t['log2ratio_501_pc1'],'dnase_2001':t['log2ratio_2001_pc1'],'adult_rna_501':v['adult_rna_log2_ratio_501bp'],'embryonic_plus_501':v['embryonic_plus_rna_log2_ratio_501bp'],'embryonic_minus_501':v['embryonic_minus_rna_log2_ratio_501bp'],'insomnia_pip':v['insomnia_conditional_abf_pip'],'adhd_pip':v['adhd_conditional_abf_pip'],'lee_allelic_fdr':x['lee_allelic_fdr'],'siraj_allelic_fdr':x['siraj_allelic_fdr']})
pth=Path(__file__).parent
with (pth/'crosswalk.tsv').open('w',newline='') as f:
 d=csv.DictWriter(f,fieldnames=out[0],delimiter='\t');d.writeheader();d.writerows(out)
(pth/'receipt.json').write_text(json.dumps({'inputs':{k:{'path':str(v),'sha256':hashlib.sha256(v.read_bytes()).hexdigest()} for k,v in inputs.items()},'source_assay_rows':len(a),'prediction_rows':len(p),'window_rows':len(w),'joined_rows':len(out)},indent=2)+'\n')
