"""Verify external constructs against hg38 and build native 270-base contexts.

Local source derivatives only. Source labels never determine matching/exclusions.
"""
from pathlib import Path
import csv,json,hashlib,collections,bisect,concurrent.futures,threading,time
import requests,numpy as np
ROOT=Path(__file__).resolve().parents[2];WORK=Path(__file__).resolve().parent;OUT=ROOT/'outputs/sleep_frailty_ai_v2';cache=WORK/'external_reference_contexts';cache.mkdir(exist_ok=True)
rows=list(csv.DictReader((WORK/'external_retallick_manifest_local.csv').open()));data=np.load(WORK/'external_retallick_local.npz');abc=np.array(list('ACGT'));rc=lambda s:s.translate(str.maketrans('ACGT','TGCA'))[::-1]
source=list(csv.DictReader((ROOT/'outputs/sleep_frailty_ai_v1/split_manifest.csv').open()));coords=collections.defaultdict(list)
for r in source:coords[r['chromosome']].append(int(r['variant'].split('-')[1]))
for c in coords:coords[c].sort()
local=threading.local()
def fetch(item):
    i,r=item;chrom='chr'+r['chromosome'];pos=int(r['source_position_build_unresolved']);start=pos-1-135;target=cache/f'{chrom}_{pos}.json'
    if target.exists():j=json.loads(target.read_text())
    else:
        if not hasattr(local,'session'):local.session=requests.Session()
        for attempt in range(3):
            try:
                q=local.session.get('https://api.genome.ucsc.edu/getData/sequence',params={'genome':'hg38','chrom':chrom,'start':start,'end':start+270},timeout=30);q.raise_for_status();j=q.json();assert len(j['dna'])==270
                j['requested_url']=q.url;j['verified_start_0based']=start;tmp=target.with_suffix('.tmp');tmp.write_text(json.dumps(j));tmp.replace(target);break
            except Exception:
                if attempt==2:raise
                time.sleep(attempt+1)
    s=j['dna'].upper()
    if len(s)!=270 or set(s)-set('ACGT'):return i,None,'unusable_native_sequence'
    a=''.join(abc[data['ref'][i]-7]);b=''.join(abc[data['alt'][i]-7]);p=int(r['changed_position']);offset=135-p
    native=s[offset:offset+199]
    if native in [a,b]:x,y=a,b
    elif s[135-(198-p):135-(198-p)+199] in [rc(a),rc(b)]:x,y=rc(a),rc(b)
    else:return i,None,'199bp_construct_flanks_do_not_match_hg38'
    # Use original source bases in the normal forward case, complement on reverse.
    if native in [a,b]:alleles={a[p],b[p]}
    else:alleles={rc(a)[198-p],rc(b)[198-p]}
    if s[135] not in alleles:return i,None,'reference_base_not_in_pair'
    alt=next(iter(alleles-{s[135]}));alternate=s[:135]+alt+s[136:]
    cc=coords[r['chromosome']];at=bisect.bisect_left(cc,pos);near=[abs(pos-cc[k]) for k in [at-1,at] if 0<=k<len(cc)];distance=min(near) if near else None
    result={'variant':r['variant'],'source_row':i,'chromosome':chrom,'position_hg38':pos,'genome_reference_allele':s[135],'tested_alternate_allele':alt,'reference_sequence':s,'alternate_sequence':alternate,'label':int(r['label']),'nearest_salomon_variant_bp':distance,'eligible_external_region':distance is None or distance>1000,'context_sha256':hashlib.sha256(target.read_bytes()).hexdigest()}
    return i,result,None
results=[];errors=[];excluded=collections.Counter()
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
    fut={pool.submit(fetch,x):x[0] for x in enumerate(rows)}
    for j,f in enumerate(concurrent.futures.as_completed(fut),1):
        try:
            i,r,reason=f.result()
            if r:results.append(r)
            else:excluded[reason]+=1
        except Exception as e:errors.append({'row':fut[f],'error':repr(e)})
        if j%50==0:print('Mapped external',j,'/',len(rows),'errors',len(errors),flush=True)
results.sort(key=lambda x:x['source_row']);retained=[r for r in results if r['eligible_external_region']]
excluded['within_1kb_of_any_salomon_variant']=len(results)-len(retained)
for name,rr in [('external_verified_all_local.csv',results),('external_pairs_local.csv',retained)]:
    if not rr:continue
    with (WORK/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rr[0]));w.writeheader();w.writerows(rr)
audit={'source_pairs':len(rows),'full_construct_hg38_matches':len(results),'retained_after_1kb_region_exclusion':len(retained),'detected_calls_retained':sum(r['label'] for r in retained),'excluded':dict(excluded),'retrieval_errors':errors,'build_verification':'Per-variant 199-base construct including an observed allele exactly matches public UCSC hg38; no coordinate-build assumption alone','context_definition':'Native 270-base hg38 window centered at index135, with single alternate allele; 199-base measured construct verified within window','limitations':['Extra native context outside the 199-base experimental construct was not tested by that experiment','1kb region exclusion is not proof against homologous sequence, training-data exposure or all regional dependencies','Source biological state and assay selection differ','External labels have not entered model fitting or selection'],'sequence_level_derivatives':'Local research only; do not redistribute'}
(OUT/'external_mapping_audit.json').write_text(json.dumps(audit,indent=2)+'\n');print(json.dumps(audit,indent=2),flush=True)
if errors:raise RuntimeError('Some public sequence queries failed; retry resumably before evaluation')
