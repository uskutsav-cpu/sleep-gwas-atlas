import csv,math,ast,json
from pathlib import Path
root=Path(__file__).resolve().parents[1];ns={'math':math};tree=ast.parse((root/'work/stage3_smim32_coloc_abf.py').read_text());exec(compile(ast.Module(body=[x for x in tree.body if isinstance(x,ast.FunctionDef)],type_ignores=[]),'coloc','exec'),ns)
a={x['variant_GRCh38_ref_alt'] for x in csv.DictReader((root/'outputs/stage5_smad5_muscle_eur_ld_variant_order.csv').open())}
def cv(x):
 c,p,r,t=x.split(':');return 'chr'+c+':'+p+':'+r+'>'+t
out=[]
for trait,sd in [('Mak_FI',1.1004195728364333),('Mak_FP',.1232461359076205),('FinnGen_insomnia',.2)]:
 rows=list(csv.DictReader((root/f'outputs/stage6_tenk10k_peak_pal_freq_{trait}_aligned.csv').open()));mask=[cv(x['variant_GRCh38']) in a for x in rows];v={'trait':trait,'total':len(rows),'in_ld':sum(mask)}
 for source,bs,ss,prior in [('qtl','qtl_beta_alt','qtl_se',.15),('gwas','gwas_beta_alt','gwas_se',sd)]:
  w=[ns['bf'](float(x[bs]),float(x[ss]),prior) for x in rows];v[source+'_single_signal_abf_support_fraction']=math.exp(ns['logsum']([b for b,m in zip(w,mask) if m])-ns['logsum'](w));i=max(range(len(w)),key=lambda i:w[i]);v[source+'_top_variant']=rows[i]['variant_GRCh38'];v[source+'_top_in_ld']=mask[i]
 out.append(v)
(root/'outputs/stage6_tenk10k_ld_support_coverage.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
