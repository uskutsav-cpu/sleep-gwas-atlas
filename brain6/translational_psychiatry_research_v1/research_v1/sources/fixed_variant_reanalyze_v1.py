#!/usr/bin/env python3
"""Process the frozen acquisition, correcting chr schema alias only.

The frozen protocol/acquisition script are preserved. No variant, ALT, endpoint,
threshold, missing rule or hypothesis changed. Outputs have corrected suffixes.
"""
import json,math
import fixed_variant_finngen as f

def main():
 m=json.loads((f.DEST/'FIXED_VARIANT_FINGEN_MANIFEST.json').read_text());native=[];by={};qcs=[]
 needed={r:{a for _,rr,s,o in f.SLOTS if rr==r for a in [s,o]} for r in f.VARIANTS}
 for v in m['variants']:
  file=f.RAW/(v['rsid']+'.json');d=json.loads(file.read_text()) if file.exists() else {}
  vv=d.get('variant',{})
  identity=all(str(vv.get(k))==str(v[x]) for k,x in [('chr','chrom'),('pos','pos38'),('ref','ref'),('alt','alt')])
  annotation=vv.get('annotation',{});rsids=annotation.get('rsids') if isinstance(annotation,dict) else None
  qc=dict(region=v['region'],rsid=v['rsid'],key=v['key'],identity_pass=identity,native_varid=vv.get('varid'),annotation_rsid=rsids,result_count=len(d.get('results',[])),source_sha256=f.digest(file) if file.exists() else '',coverage='EXACT_POINT_ONLY; dense_region_completeness_UNKNOWN')
  qcs.append(qc)
  for code in sorted(needed[v['region']]):
   rr=[x for x in d.get('results',[]) if x.get('phenocode')==code];z=rr[0] if len(rr)==1 and identity else {}
   reason='ESTIMATED' if z else ('IDENTITY_FAILED' if d and not identity else 'MISSING_EXACT_POINT_OR_PHENOTYPE')
   required=all(f.finite(z.get(k)) for k in ['pval','beta','sebeta','maf'])
   required=required and 0<z.get('pval',-1)<=1 and z.get('sebeta',-1)>0 and 0<=z.get('maf',-1)<=1
   if z and not required:reason='INVALID_OR_MISSING_REQUIRED_NATIVE_FIELD'
   row=dict(region=v['region'],rsid=v['rsid'],key=v['key'],phenocode=code,status=reason,source_sha256=qc['source_sha256'],identity_pass=identity,effect_allele=v['alt'],other_allele=v['ref'])
   for k in ['phenostring','pval','beta','sebeta','maf','maf_case','maf_control','n_case','n_control','n_sample','mlogp']:row[k]=z.get(k)
   if reason=='ESTIMATED':
    row.update(odds_ratio=math.exp(z['beta']),or_ci95_lower=math.exp(z['beta']-1.959963984540054*z['sebeta']),or_ci95_upper=math.exp(z['beta']+1.959963984540054*z['sebeta']),wald_p_diagnostic=math.erfc(abs(z['beta']/z['sebeta'])/math.sqrt(2)))
    nc,nn=z.get('n_case'),z.get('n_control');row['neff_count_formula']=4*nc*nn/(nc+nn) if f.finite(nc) and f.finite(nn) and nc>0 and nn>0 else None
   else:row.update(odds_ratio=None,or_ci95_lower=None,or_ci95_upper=None,wald_p_diagnostic=None,neff_count_formula=None)
   row['sample_size_semantics']='Native endpoint metadata/count fields; not asserted per-variant observed N';native.append(row);by[v['region'],code]=row
 slots=[]
 for s,r,a,b in f.SLOTS:
  x,y=by[r,a],by[r,b];ok=x['status']==y['status']=='ESTIMATED';pc=max(x['pval'],y['pval']) if ok else None
  slots.append(dict(slot_id=s,region=r,rsid=f.VARIANTS[r]['rsid'],sleep_code=a,outcome_code=b,sleep_status=x['status'],outcome_status=y['status'],sleep_native_p=x['pval'],outcome_native_p=y['pval'],conjunction_p=pc,bonferroni_p=min(1,4*pc) if ok else None,alpha_per_slot=f.ALPHA,decision=('BOTH_ASSOCIATIONS_AT_FIXED_VARIANT' if pc<=f.ALPHA else 'DID_NOT_MEET_FAMILY_THRESHOLD') if ok else 'NOT_ESTIMATED',classification='PHENOTYPE_SENSITIVITY',cohort_independence='UNVERIFIED',same_causal_variant='NOT_TESTED',direction_relative_discovery='NOT_A_FROZEN_ENDPOINT'))
 f.tsv(f.DEST/'fixed_variant_native_fields_corrected.tsv',native);f.tsv(f.DEST/'fixed_variant_results_corrected.tsv',slots);f.tsv(f.DEST/'fixed_variant_qc_corrected.tsv',qcs)
 print(json.dumps(slots,indent=2))
if __name__=='__main__':main()
