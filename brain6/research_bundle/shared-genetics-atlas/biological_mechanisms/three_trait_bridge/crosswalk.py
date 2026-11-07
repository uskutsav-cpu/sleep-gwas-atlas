import csv,hashlib,json
from pathlib import Path
base=Path(__file__).resolve().parents[2]; dst=Path(__file__).parent
inputs={'reporter':base/'biological_mechanisms/cross_assay_prediction/lee_27_crosswalk.tsv','scores':base/'variant_context/brain6_fi_ffs_screen/three_score.tsv','symptom':Path('/Volumes/Extreme SSD/brain6-three-trait-candidates-results-20261028.tsv'),'pred':base.parent/'brain6-alphagenome-pilot/results/jbrowse_expanded/expanded_variant_evidence.tsv'}
def read(p):
 with p.open(newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
rep=read(inputs['reporter']); score=read(inputs['scores']); three=read(inputs['symptom']); pred=read(inputs['pred'])
def idx(a,key):
 d={x[key]:x for x in a}; assert len(d)==len(a);return d
r=idx(rep,'rsid'); f=idx(score,'variant'); g=idx(three,'rsid'); p=idx(pred,'snp'); assert len(r)==27 and len(f)==len(g)==len(p)==425
selected=set(r)|{'rs62264764','rs1476535','rs4916723','rs2431108','rs77960'}
out=[]
for rs in sorted(selected):
 a=p[rs];b=f[rs];c=g[rs];assert a['candidate_locus_id']==b['region']==c['region'] and a['atlas_variant']==c['grch38'];x=r.get(rs,{})
 out.append({'rsid':rs,'region':b['region'],'grch38':c['grch38'],'selected_reporter_27':rs in r,'lee_class':x.get('lee_allelic_class',''),'lee_allelic_fdr':x.get('lee_allelic_fdr',''),'adult_dnase_501':a['dnase_log2_ratio_501bp'],'insomnia_pip':a['insomnia_conditional_abf_pip'],'adhd_pip':a['adhd_conditional_abf_pip'],'p_fi':b['p_FI'],'p_ffs':b['p_FFS'],'p_ukb_hfrs':b['p_HFRS'],'p_mvp_insomnia':c['mvp_p'],'p_finn_hfrs':c['hfrs_p'],'p_childhood_symptoms':c['symptom_p'],'three_trait_maxp_adjusted_425':c['joint_bonferroni_425']})
with (dst/'crosswalk.tsv').open('w',newline='') as fp:
 w=csv.DictWriter(fp,fieldnames=out[0],delimiter='\t');w.writeheader();w.writerows(out)
(dst/'receipt.json').write_text(json.dumps({'input_sha256':{k:hashlib.sha256(v.read_bytes()).hexdigest() for k,v in inputs.items()},'selected_variants':len(out),'source_reporter_rows':len(rep),'scores':len(score),'three_trait':len(three),'predictions':len(pred)},indent=2)+'\n')
for x in out:
 if x['rsid'] in {'rs62264764','rs1476535','rs4916723','rs2431108','rs77960','rs325485','rs6874138','rs62369151'}:print(x)
