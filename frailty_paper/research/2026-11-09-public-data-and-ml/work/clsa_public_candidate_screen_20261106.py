"""Cache exact public CLSA lookups and audit three physical-component outcomes.

No individual records or full-score prediction are analyzed. Outcomes were selected
from their labels before reading the 33 candidate association results. The first
MEIS1 lookup was inspected during API discovery; this is explicitly exploratory.
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import csv, json, hashlib, subprocess, datetime

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/clsa_public_followup_20261106'
CACHE=OUT/'raw_variant_responses'
CACHE.mkdir(parents=True,exist_ok=True)
BASE='https://clsa-pheweb.cerc-genomic-medicine.ca/api/v1'
TARGETS=['GS_EXAM_AVG_COM','GS_EXAM_MAX_COM','WLK_TIME_COM']

def readcsv(path):
    with path.open(newline='') as f:return list(csv.DictReader(f))

def writecsv(path,rows,fields=None):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields or list(rows[0])); w.writeheader();w.writerows(rows)

catalog_path=ROOT/'work/clsa_phenotypes_20261106.json'
catalog=json.loads(catalog_path.read_text())
euro=[r for r in catalog if r['stratification']=={'ancestry':'european','sex':'combined'}]
writecsv(OUT/'target_phenotypes.csv',[{k:r[k] for k in ['phenocode','phenostring','num_samples','num_cases','num_controls','category']} for r in euro if r['phenocode'] in TARGETS])
rows=readcsv(ROOT/'outputs/ncds_score_variant_request.csv')
variants={}
for r in rows:
    if r['unique_reciprocal_position']!='True':raise ValueError(r)
    k=r['variant']; pair=tuple(sorted([r['counted_allele_38_candidate'],r['other_allele_38_candidate']]))
    if k in variants:
        assert variants[k]['pair']==pair and variants[k]['pos38']==r['pos38']
        variants[k]['score_membership'].append(r['score'])
    else:variants[k]={'variant':k,'chrom38':r['chr38'],'pos38':r['pos38'],'pair':pair,'score_membership':[r['score']]}

def query(r):
    errors=[]
    for a,b in [r['pair'],r['pair'][::-1]]:
        vid=f"{r['chrom38']}-{r['pos38']}-{a}-{b}"
        path=CACHE/f'{r["variant"]}_{vid}.json'
        url=f'{BASE}/variant/{vid}/european.combined'
        try:
            if path.exists():raw=path.read_bytes()
            else:
                # curl uses the host's configured certificate trust store.
                raw=subprocess.run(['curl','--fail','--location','--silent','--show-error','--max-time','40',url],check=True,capture_output=True).stdout
                path.write_bytes(raw)
            d=json.loads(raw)
            if not isinstance(d,dict) or 'phenos' not in d:
                errors.append({'url':url,'error':str(d)[:180]});continue
            assert str(d['chrom'])==r['chrom38'] and str(d['pos'])==r['pos38']
            assert {d['ref'],d['alt']}==set(r['pair'])
            # Browser may list multiple rsIDs; require the original identifier.
            ids=str(d.get('rsids','')).replace(';',',').replace(' ',',').split(',')
            assert r['variant'] in ids,(r['variant'],ids)
            return {'variant':r['variant'],'status':'exact_coordinate_allele_pair_and_rsid','url':url,'response_sha256':hashlib.sha256(raw).hexdigest(),'n_phenos':len(d['phenos']),'data':d,'errors':errors}
        except Exception as exc:errors.append({'url':url,'error':str(exc)[:200]})
    return {'variant':r['variant'],'status':'unresolved','errors':errors,'data':None}

with ThreadPoolExecutor(max_workers=4) as pool:results=list(pool.map(query,variants.values()))
receipt={r['variant']:{k:v for k,v in r.items() if k!='data'} for r in results}
matrix=[]
for r in results:
    v=variants[r['variant']];d=r['data']
    for ph in TARGETS:
        matches=[] if d is None else [x for x in d['phenos'] if x['phenocode']==ph]
        assert len(matches)<=1
        x=matches[0] if matches else {}
        matrix.append({'variant':r['variant'],'chr38':v['chrom38'],'pos38':v['pos38'],'ref':d.get('ref','') if d else '', 'alt_effect_allele_browser':d.get('alt','') if d else '', 'phenocode':ph,'phenostring':x.get('phenostring',''),'num_samples_catalog':x.get('num_samples',''),'n_samples_variant':x.get('n_samples',''),'beta_alt_rounded':x.get('beta',''),'se_alt_rounded':x.get('sebeta',''),'p_value_rounded':x.get('pval',''),'alt_frequency_rounded':x.get('af',''),'imputation_quality_rounded':x.get('imp_quality',''),'test':x.get('test',''),'status':'matched' if x else 'unavailable','source_url':r.get('url','')})

# Holm correction across the entire planned 33 x 3 family, retaining unavailable
# tests in the family size. Source browser precision bounds all interpretations.
ordered=sorted([r for r in matrix if r['status']=='matched'],key=lambda r:r['p_value_rounded'])
previous=0
for i,r in enumerate(ordered):
    previous=max(previous,min(1,(len(matrix)-i)*r['p_value_rounded']))
    r['p_holm_99_browser_rounded']=previous
for r in matrix:r.setdefault('p_holm_99_browser_rounded','')
writecsv(OUT/'candidate_component_associations.csv',matrix)
aligned=[]
for name in ['unvalidated_15_variant_score_prototype.csv','unvalidated_22_variant_score_prototype.csv','unvalidated_22_shared_direction_score.csv']:
    scores=readcsv(ROOT/'outputs'/name)
    for s in scores:
        variant=s.get('variant',s.get('rsid'))
        # Use per-model reciprocal-coordinate allele mapping, not an rsID-only flip.
        # The actual manifest uses 22_count for the equal-count convention.
        mapped=[x for x in rows if x['variant']==variant and ((len(scores)==15 and x['score']=='15_weighted') or (len(scores)==22 and 'prototype' in name and x['score']=='22_weighted') or (len(scores)==22 and 'shared_direction' in name and x['score'] not in ['15_weighted','22_weighted']))]
        assert len(mapped)==1,(name,variant,mapped)
        m=mapped[0]
        for x in matrix:
            if x['variant']!=variant or x['status']!='matched':continue
            counted=m['counted_allele_38_candidate']; other=m['other_allele_38_candidate']
            assert {counted,other}=={x['ref'],x['alt_effect_allele_browser']}
            sign=1 if counted==x['alt_effect_allele_browser'] else -1
            aligned.append({'score_file':name,'variant':variant,'phenocode':x['phenocode'],'counted_allele_hg38':counted,'other_allele_hg38':other,'beta_counted_rounded':sign*x['beta_alt_rounded'],'se_rounded':x['se_alt_rounded'],'p_value_rounded':x['p_value_rounded'],'p_holm_99_browser_rounded':x['p_holm_99_browser_rounded'],'contrast_warning':'Pair matches target; source multiallelic comparison group and source FFS transformation remain unresolved; no predictive validation'})
writecsv(OUT/'score_allele_aligned_component_associations.csv',aligned)
summary={'checked_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'variants':len(variants),'targets':TARGETS,'planned_tests':len(matrix),'matched_tests':len(ordered),'matched_variants':sum(r['data'] is not None for r in results),'nominal_p_below_05':sum(r['p_value_rounded']<.05 for r in ordered),'holm_below_05':sum(r['p_holm_99_browser_rounded']<.05 for r in ordered),'catalog_rows':len(catalog),'catalog_unique_phenotypes':len({r['phenocode'] for r in catalog}),'catalog_frailty_or_insomnia_labels':[r['phenocode'] for r in catalog if any(z in (r['phenostring']+' '+r['phenocode']).lower() for z in ['frail','insomn'])],'catalog_sha256':hashlib.sha256(catalog_path.read_bytes()).hexdigest(),'limitations':'Exploratory rounded public marginal association lookups, not a participant-level score, a complete frailty phenotype, causal inference, or proven discovery/relative independence. Genome build and effect ALT convention supported by public frontend; original source contrasts remain subject to audit.'}
(OUT/'receipt.json').write_text(json.dumps({'summary':summary,'lookups':receipt},indent=2)+'\n')
print(json.dumps(summary,indent=2))
for r in ordered[:8]:print(r['variant'],r['phenocode'],r['beta_alt_rounded'],r['p_value_rounded'],r['p_holm_99_browser_rounded'])
