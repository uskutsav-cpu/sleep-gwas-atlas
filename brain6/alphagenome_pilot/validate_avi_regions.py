#!/usr/bin/env python3
"""Cross-check one prespecified high-AVI variant per region with Broad's lookup.

The Broad table and mirror derive from the same AlphaGenome data: agreement is
an independent data-access check, not independent biological validation.
"""
import csv,hashlib,json,time,urllib.parse,urllib.request,ssl
from pathlib import Path
from datetime import datetime,timezone
P=Path(__file__).resolve().parent/'results'/'avi_mirror'; INPUT=P/'variant_evidence_avi.tsv'; OUT=P/'broad_crosscheck.tsv'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
with INPUT.open() as f: rows=list(csv.DictReader(f,delimiter='\t'))
groups={}
for r in rows:
    if r['avi_status']!='SCORED_MIRROR_EXPLORATORY':raise ValueError('Unscored row')
    groups.setdefault(r['candidate_locus_id'],[]).append(r)
assert len(groups)==6
chosen=[sorted(rs,key=lambda r:(-int(r['avi_raw_e5']),r['atlas_variant']))[0] for _,rs in sorted(groups.items())]
ctx=ssl.create_default_context(cafile='/etc/ssl/cert.pem'); out=[]
for i,r in enumerate(chosen):
    if i:time.sleep(12) # At most a handful of requests per minute.
    params={'hg':'38','chrom':r['chromosome'],'pos':r['pos_grch38_1based'],'ref':r['ref_grch38'],'alt':r['alt_grch38']}
    url='https://score-lookup-xwkwwwxdwq-uc.a.run.app/?'+urllib.parse.urlencode(params)
    result={'candidate_locus_id':r['candidate_locus_id'],'snp':r['snp'],'atlas_variant':r['atlas_variant'],'lookup_url':url,'mirror_raw':r['avi_raw_score_decimal'],'mirror_phred_low':r['avi_phred_lower'],'mirror_phred_high':r['avi_phred_upper'],'broad_raw':'','broad_phred':'','status':'','detail':''}
    try:
        with urllib.request.urlopen(url,context=ctx,timeout=25) as response:data=json.load(response)
        v=data.get('variant',{}); avi=data.get('scores',{}).get('alphagenome_avi'); errors=data.get('errors',[])
        if (v.get('hg'),v.get('chrom'),str(v.get('pos')),v.get('ref'),v.get('alt'))!=('38',params['chrom'],params['pos'],params['ref'],params['alt']):raise ValueError('Returned variant identity mismatch')
        if any(e.get('table')=='AlphaGenome AVI' for e in errors):result.update(status='SERVICE_ERROR',detail='AVI table lookup error')
        elif not avi:result.update(status='NOT_RETURNED',detail='No AVI score in response')
        else:
            raw=float(avi['rawScore']);phred=float(avi['phredScore']);result['broad_raw']=raw;result['broad_phred']=phred
            agree_raw=abs(raw-float(r['avi_raw_score_decimal']))<=0.000501
            agree_phred=float(r['avi_phred_lower'])-0.00501<=phred<=float(r['avi_phred_upper'])+0.00501
            result['status']='AGREES_TO_ROUNDING' if agree_raw and agree_phred else 'MISMATCH'
    except Exception as e:result.update(status='QUERY_ERROR',detail=type(e).__name__)
    out.append(result)
    with OUT.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(out[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(out)
    print(f"{i+1}/6 {result['snp']} {result['status']}",flush=True)
receipt={'source_evidence_sha256':sha(INPUT),'output_sha256':sha(OUT),'selection':'highest raw AVI per region, ties by GRCh38 variant string; chosen before lookup','source':'https://github.com/broadinstitute/SpliceAI-lookup/tree/main/google_cloud_run_services/score_lookup','timestamp_utc':datetime.now(timezone.utc).isoformat(),'n':len(out),'statuses':{s:sum(r['status']==s for r in out) for s in sorted({r['status'] for r in out})},'tolerance':'raw 0.000501 and PHRED 0.00501 plus mirror knot interval; Broad outputs rounded','scope':'Independent lookup and formatting verification from the same source data; not biological validation; not verification of remaining variants.'}
(P/'broad_crosscheck_provenance.json').write_text(json.dumps(receipt,indent=2)+'\n')
