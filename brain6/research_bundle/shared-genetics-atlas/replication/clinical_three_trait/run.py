import csv,hashlib,json,math
from pathlib import Path
here=Path(__file__).resolve().parent; adhd=here.parent/'finngen_adhd_raw_lookup/indexed_results.tsv'; pair=Path('/Volumes/Extreme SSD/brain6-mvp-finngen-hfrs-results-20261028.tsv'); mvp=here.parent/'mvp_finngen_clinical/results.tsv'
def read(p):
 with p.open(newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
a=read(adhd); b=read(pair); m={x['variant']:x for x in read(mvp)}; assert len(a)==7 and len(b)==425; by={x['rsid']:x for x in b}; assert len(by)==425
out=[]
for x in a:
 y=by.get(x['variant']); row={'variant':x['variant'],'clinical_status':x['status'],'candidate_status':'not_in_425' if y is None else 'matched_425','grch38':y['grch38'] if y else '', 'mvp_insomnia_p':m.get(x['variant'],{}).get('mvp_insomnia_p',''),'finngen_adhd_p':x['source_p'],'finngen_hfrs_p':'','adhd_beta_alt':'','insomnia_beta_alt':'','hfrs_beta_alt':'','max_p':'','max_p_times_7':'','status':'not_estimable'}
 if y and x['status']=='MATCHED' and y['mvp_status']=='matched' and y['hfrs_status']=='matched':
  key=y['grch38'];prefix,alleles=key.split(':',2)[0:2],key.split(':',2)[2];ref,alt=alleles.split('>');assert prefix[0]=='chr'+x['chr'] and int(prefix[1])==int(x['position38']); assert {ref,alt}=={x['effect_allele'],x['other_allele']}
  ps=[float(y['mvp_p']),float(x['source_p']),float(y['hfrs_p'])]; assert all(math.isfinite(z) and 0<=z<=1 for z in ps)
  eff=float(x['beta_aligned'])*(1 if x['effect_allele']==alt else -1)
  row.update({'mvp_insomnia_p':ps[0],'finngen_hfrs_p':ps[2],'adhd_beta_alt':eff,'insomnia_beta_alt':y['mvp_beta_alt'],'hfrs_beta_alt':y['hfrs_beta_alt'],'max_p':max(ps),'max_p_times_7':min(1,7*max(ps)),'status':'matched'})
 out.append(row)
assert len(out)==7
with (here/'results.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=out[0],delimiter='\t');w.writeheader();w.writerows(out)
(here/'receipt.json').write_text(json.dumps({'inputs_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in (adhd,pair,mvp)},'matched':sum(x['status']=='matched' for x in out),'total':len(out)},indent=2)+'\n')
for x in out:print(x)
