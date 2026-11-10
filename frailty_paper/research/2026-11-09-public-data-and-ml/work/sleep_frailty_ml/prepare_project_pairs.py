"""Fetch public hg38 reference contexts for the existing 33 candidate contrasts."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import json,csv,subprocess,urllib.parse,hashlib

ROOT=Path(__file__).resolve().parents[2];WORK=Path(__file__).resolve().parent
OUT=ROOT/'outputs/sleep_frailty_ai_v1';CACHE=WORK/'reference_contexts';CACHE.mkdir(exist_ok=True)
with (ROOT/'outputs/ncds_score_variant_request.csv').open() as f:rows=list(csv.DictReader(f))
unique={}
for r in rows:
    assert r['unique_reciprocal_position']=='True'
    pair=sorted([r['counted_allele_38_candidate'],r['other_allele_38_candidate']])
    if r['variant'] in unique:assert unique[r['variant']]['allele_pair']==pair
    else:unique[r['variant']]={'variant':r['variant'],'chromosome':r['chr38'],'position':int(r['pos38']),'allele_pair':pair}

def fetch(r):
    start=r['position']-1-134;end=start+270
    url='https://api.genome.ucsc.edu/getData/sequence?'+urllib.parse.urlencode({'genome':'hg38','chrom':'chr'+r['chromosome'],'start':start,'end':end})
    file=CACHE/f'{r["variant"]}.json'
    try:
        if file.exists():raw=file.read_bytes()
        else:raw=subprocess.run(['curl','-fLsS','--max-time','35',url],check=True,capture_output=True).stdout;file.write_bytes(raw)
        d=json.loads(raw);s=d['dna'].upper();assert len(s)==270 and not set(s)-set('ACGT')
        assert d['start']==start and d['end']==end and d['chrom']=='chr'+r['chromosome']
        assert s[134] in r['allele_pair'],(s[134],r['allele_pair'])
        other=next(a for a in r['allele_pair'] if a!=s[134]);assert len(other)==1
        alt=s[:134]+other+s[135:]
        return {'variant':r['variant'],'chromosome':r['chromosome'],'position_hg38':r['position'],'genome_reference_allele':s[134],'tested_alternate_allele':other,'reference_sequence':s,'alternate_sequence':alt,'source_url':url,'source_sha256':hashlib.sha256(raw).hexdigest(),'status':'reference_base_and_exact_target_pair_verified'}
    except Exception as e:return {'variant':r['variant'],'status':'unresolved','error':str(e),'source_url':url}

with ThreadPoolExecutor(max_workers=4) as pool:results=list(pool.map(fetch,unique.values()))
valid=[r for r in results if r['status'].startswith('reference_base')]
with (WORK/'project_sequence_pairs.csv').open('w',newline='') as f:
    writer=csv.DictWriter(f,fieldnames=list(valid[0]));writer.writeheader();writer.writerows(valid)
manifest=[{k:v for k,v in r.items() if 'sequence' not in k} for r in results]
(OUT/'project_reference_context_receipt.json').write_text(json.dumps(manifest,indent=2)+'\n')
print('Verified project contexts:',len(valid),'of',len(unique),flush=True)
for r in results:
    if r['status']=='unresolved':print(r,flush=True)
