import csv,hashlib,json,math
from pathlib import Path
from statistics import NormalDist
p=Path(__file__).resolve().parent;root=p.parents[1]
paths=[root/'mvp_six_leads/six_leads.tsv',root/'mvp_six_leads/receipt.json',root/'six_lead_joint/six_leads.tsv']
def read(p):
 with p.open(newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
m=read(paths[0]);f={x['variant']:x for x in read(paths[2])};a={x['variant']:x for x in json.loads(paths[1].read_text())['audit']}
assert len(m)==6 and len(f)==6
out=[]
for x in m:
 y=f[x['variant']];assert x['chr']==y['chr'] and x['pos38']==y['position38'] and x['effect_allele']==y['effect_allele']
 for method,rho in [(method,rho) for method in ('ci_derived','p_derived') for rho in (-0.2,0,0.2,0.5)]:
  r={'variant':x['variant'],'mvp_se_method':method,'rho_assumed':rho,'status':'missing','mvp_beta':'','mvp_se_approx':'','mvp_p_original':'','mvp_p_approx':'','finngen_beta':'','finngen_se':'','difference_mvp_minus_finngen':'','difference_se':'','ci_low':'','ci_high':'','p_difference':'','p_times_six':''}
  if y['status']=='MATCHED' and x['status']=='MATCHED':
   q=a[x['variant']]['mvp_exact_rows'][0];lo=float(q['ci_lower']);hi=float(q['ci_upper']);assert 0<lo<hi
   b=float(x['mvp_logor_aligned']);se_ci=(math.log(hi)-math.log(lo))/(2*NormalDist().inv_cdf(.975));se_p=abs(b)/NormalDist().inv_cdf(1-float(q['p_value'])/2);se=se_ci if method=='ci_derived' else se_p;fb=float(y['sleep_beta_aligned']);fs=float(y['sleep_se']);dif=b-fb;dse=math.sqrt(se**2+fs**2-2*rho*se*fs);z=dif/dse;pv=math.erfc(abs(z)/math.sqrt(2));
   r.update(status='matched',mvp_beta=b,mvp_se_approx=se,mvp_p_original=q['p_value'],mvp_p_approx=math.erfc(abs(b/se)/math.sqrt(2)),finngen_beta=fb,finngen_se=fs,difference_mvp_minus_finngen=dif,difference_se=dse,ci_low=dif-1.95996398454*dse,ci_high=dif+1.95996398454*dse,p_difference=pv,p_times_six=min(1,pv*6))
  out.append(r)
with (p/'effects.tsv').open('w',newline='') as h:
 w=csv.DictWriter(h,fieldnames=out[0],delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(out)
(p/'receipt.json').write_text(json.dumps({'inputs_sha256':{str(t.relative_to(root)):hashlib.sha256(t.read_bytes()).hexdigest() for t in paths},'n_markers':6,'matched_markers':sum(r['status']=='matched' and r['rho_assumed']==0 and r['mvp_se_method']=='ci_derived' for r in out),'assumed_rho':[-.2,0,.2,.5],'mvp_se':'CI and p derived approximations; source test and true SE not verified','note':'Postselected; overlap unknown; different phenotypes and populations'},indent=2)+'\n')
for r in out:
 if r['rho_assumed']==0:print(r['mvp_se_method'],r['variant'],r['status'],r['difference_mvp_minus_finngen'],r['p_difference'],r['p_times_six'],r['mvp_p_approx'])
