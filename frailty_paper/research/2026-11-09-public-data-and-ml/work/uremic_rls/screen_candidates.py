#!/usr/bin/env python3
"""Extract exact project allele contrasts from publicly released ESRD/RLS GWAS."""
import csv, hashlib, json, math, os
from collections import defaultdict

base='work/uremic_rls/'
manifest=json.load(open(base+'source_manifest.json'))
models=defaultdict(list)
with open('outputs/ncds_score_variant_request.csv',newline='') as f:
    for r in csv.DictReader(f):
        key=(str(r['chr37']),int(r['pos37']))
        entry=(r['variant'],r['counted_allele_37'],r['other_allele_37'],r['score'])
        if entry not in models[key]: models[key].append(entry)
names=['GER1','GER2','GRK']; out=[]; stats=[]
for study in names:
    meta=next(x for x in manifest['files'] if x['name'].startswith(study+'_GWAS'))
    path=base+study+'.txt'
    if os.path.getsize(path)!=meta['size']: raise RuntimeError(f'{study}: unexpected size')
    h=hashlib.md5()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    if 'md5:'+h.hexdigest()!=meta['checksum']: raise RuntimeError(f'{study}: checksum mismatch')
    seen=0; matched=0; invalid=0
    with open(path,newline='') as f:
        reader=csv.DictReader(f,delimiter='\t')
        required={'chromosome','base_pair_location','effect_allele','other_allele','beta','standard_error','effect_allele_frequency','p_value'}
        if not required.issubset(reader.fieldnames): raise RuntimeError('Unexpected schema')
        for r in reader:
            seen+=1; key=(r['chromosome'].removeprefix('chr'),int(r['base_pair_location']))
            if key not in models:continue
            for variant,counted,other,score in models[key]:
                a1=r['effect_allele'].upper();a2=r['other_allele'].upper()
                if {a1,a2}!={counted,other}: continue
                try:
                    b=float(r['beta']); se=float(r['standard_error']); p=float(r['p_value']); eaf=float(r['effect_allele_frequency'])
                except ValueError:
                    invalid+=1; continue
                if not all(math.isfinite(x) for x in [b,se,p,eaf]) or not (se>0 and 0<=p<=1 and 0<=eaf<=1):
                    invalid+=1; continue
                matched+=1
                out.append({'study':study,'variant':variant,'score':score,'chr37':key[0],'pos37':key[1],
                            'counted_allele':counted,'other_allele':other,'source_effect_allele':a1,
                            'beta_counted':b if a1==counted else -b,'se':se,'p':p,
                            'counted_frequency':eaf if a1==counted else 1-eaf})
    stats.append({'study':study,'input_md5':h.hexdigest(),'rows_scanned':seen,'matched_score_rows':matched,'invalid_matched_score_rows':invalid})
with open('outputs/stage6_uremic_rls_candidate_lookup.csv','w',newline='') as f:
    writer=csv.DictWriter(f,fieldnames=['study','variant','score','chr37','pos37','counted_allele','other_allele','source_effect_allele','beta_counted','se','p','counted_frequency']);writer.writeheader();writer.writerows(out)
with open('outputs/stage6_uremic_rls_candidate_receipt.json','w') as f:
    json.dump({'files':stats,'result_rows':len(out),'source':'https://doi.org/10.5281/zenodo.21234712',
        'limitations':'GRCh37 position and stated allele pair match only; effect flip assumes same biallelic contrast. Absence is not a null result. Source releases are individual-study summary statistics, not individual genotype data; selection, multiplicity, cohort overlap and disease-specific sampling matter.'},f,indent=2)
print(json.dumps({'files':stats,'rows':len(out)},indent=2))
