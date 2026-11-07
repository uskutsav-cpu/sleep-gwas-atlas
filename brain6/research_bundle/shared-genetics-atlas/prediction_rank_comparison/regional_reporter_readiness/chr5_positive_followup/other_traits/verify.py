import hashlib,json,math,pathlib
p=pathlib.Path(__file__).resolve().parent;r=json.loads((p/'receipt.json').read_text());assert all(v['status']=='MATCHED' and len(v['exact_matches'])==1 for v in r.values())
m=r['mvp']['exact_matches'][0];h=r['hfrs']['exact_matches'][0];assert (m['effect_allele'],m['other_allele'])==('G','A') and (h['effect_allele'],h['other_allele'])==('A','G')
beta_a=-math.log(float(m['odds_ratio']));pv_m=float(m['p_value']);pv_h=float(h['p_value']);sym=json.loads((p.parent/'adhd_symptom/receipt.json').read_text())['matched'][0];pv_a=float(sym['p_value'])
out=dict(status='PASS',mvp_beta_A_from_rounded_or=beta_a,hfrs_beta_A=float(h['beta']),max_p_mvp_hfrs=max(pv_m,pv_h),max_p_mvp_hfrs_symptom=max(pv_m,pv_h,pv_a),note='post-selected raw marginal maxima, not multiplicity-corrected results',input_sha256=hashlib.sha256((p/'receipt.json').read_bytes()).hexdigest());(p/'verification.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
