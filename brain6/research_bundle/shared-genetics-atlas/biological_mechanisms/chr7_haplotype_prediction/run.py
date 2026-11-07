import csv,hashlib,json
from pathlib import Path
base=Path(__file__).resolve().parents[2]; here=Path(__file__).parent
paths={'coloc':base/'frailty_coloc/conditional_variant_posteriors.tsv','pred':base.parent/'brain6-alphagenome-pilot/results/jbrowse_expanded/expanded_variant_evidence.tsv','window':base/'alphagenome_window_sensitivity/window_scores.tsv','frailty':base/'variant_context/brain6_fi_ffs_screen/three_score.tsv'}
def read(f):
 with f.open(newline='') as h:return list(csv.DictReader(h,delimiter='\t'))
a,b,c,d=(read(paths[k]) for k in ('coloc','pred','window','frailty')); ids=['rs10953765','rs1456031','rs1476535'];out=[]
for rs in ids:
 loc=[x for x in a if x['rsid']==rs and x['chr']=='7'];pv=[x for x in b if x['snp']==rs];win=[x for x in c if x['snp']==rs];score=[x for x in d if x['variant']==rs]; assert len(loc)==len(pv)==len(win)==len(score)==1
 x=loc[0];y=pv[0];w=win[0];q=score[0];assert x['bp37']==y['pos_grch37_1based'] and y['atlas_variant']==w['variant'] and y['prediction_window_sha256']==w['input_sha256']; assert set((x['effect_allele'],x['other_allele']))==set((y['ref_grch38'],y['alt_grch38']))
 flip=1 if x['effect_allele']==y['alt_grch38'] else -1
 out.append({'rsid':rs,'variant_grch38':y['atlas_variant'],'coloc_conditional_h4':x['conditional_h4'],'fi_beta_alt':flip*float(x['fi_beta']),'finngen_insomnia_beta_alt':flip*float(x['r13_beta']),'fi_p':x['fi_p'],'finngen_insomnia_p':x['r13_p'],'ffs_p':q['p_FFS'],'ukb_hfrs_p':q['p_HFRS'],'dnase_101':w['log2ratio_101_pc1'],'dnase_501':w['log2ratio_501_pc1'],'dnase_2001':w['log2ratio_2001_pc1']})
with (here/'crosswalk.tsv').open('w',newline='') as f:
 wr=csv.DictWriter(f,fieldnames=out[0],delimiter='\t');wr.writeheader();wr.writerows(out)
(here/'receipt.json').write_text(json.dumps({'inputs_sha256':{k:hashlib.sha256(v.read_bytes()).hexdigest() for k,v in paths.items()},'rows':len(out),'condition':'conditional SNP probabilities given H4'},indent=2)+'\n')
for z in out:print(z)
