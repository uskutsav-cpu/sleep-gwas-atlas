import hashlib,json,math,pathlib
p=pathlib.Path(__file__).resolve().parent;r=json.loads((p/'receipt.json').read_text());assert r['status']=='MATCHED' and len(r['matched'])==1
x=r['matched'][0];assert x['rsid']=='rs12659965' and x['effect_allele']=='A' and x['other_allele']=='G';b=float(x['beta']);se=float(x['standard_error']);pv=float(x['p_value']);wald=math.erfc(abs(b/se)/math.sqrt(2));assert abs(wald-pv)<1e-10
parent=json.loads((p.parent/'receipt.json').read_text());sl=parent['source_reports']['R13']['rows'][0];assert (sl['ref'],sl['alt'])==('G','A');sp=float(sl['pval']);maxp=max(sp,pv)
out=dict(status='PASS',symptom_z=b/se,symptom_wald_p=wald,normal_95_low=b-1.96*se,normal_95_high=b+1.96*se,source_sleep_p=sp,unadjusted_max_p=maxp,warning='post-selected descriptive number, not a multiplicity-corrected result',input_sha256=hashlib.sha256((p/'receipt.json').read_bytes()).hexdigest())
(p/'verification.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
