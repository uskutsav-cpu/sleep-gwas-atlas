import json,requests,concurrent.futures
from pathlib import Path
rows=json.loads(Path('outputs/stage6_hansen_15_candidate_screen.json').read_text())['rows'];rows=[r for r in rows if r['exact_matches']==1]
def lookup(r):
 p=r['position_grch38'];url=f"https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr{r['chromosome']};start={p-1};end={p}"
 q=requests.get(url,timeout=30);q.raise_for_status();base=q.json()['dna'].upper();return {'variant':r['variant'],'reference':base,'source_ref':r['source_ref'],'match':base==r['source_ref']}
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as e:results=list(e.map(lookup,rows))
Path('outputs/stage6_hansen_15_reference_check.json').write_text(json.dumps(results,indent=2)+'\n');print('verified',len(results),'matches',sum(x['match'] for x in results))
