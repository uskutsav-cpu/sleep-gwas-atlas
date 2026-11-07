import csv,hashlib,json,math
from collections import Counter
from pathlib import Path
p=Path(__file__).resolve().parent
with (p/'eur_ld_blocks.bed').open() as f:
 rows=[line.split() for line in f if line.strip()]
assert rows[0]==['chr','start','stop'];blocks=[(c,int(s),int(e)) for c,s,e in rows[1:]]
assert len(blocks)==1703 and all(s<e for c,s,e in blocks)
assert len(set(blocks))==len(blocks) and len([x for x in blocks if x==('chr3',47727212,49316972)])==1
with (p.parent/'song_local_reported.tsv').open() as f:source=next(csv.DictReader(f,delimiter='\t'))
assert (int(source['Start']),int(source['End']))==(47727212,49316972)
nominal=float(source['P']); assert nominal==9.45e-6
with (p/'lava_gate_ledger.tsv').open() as f:ledger=list(csv.DictReader(f,delimiter='\t'))
assert len(ledger)==12 and all(x['frozen_gate_pass']=='False' for x in ledger)
for x in ledger:
 assert int(x['process_failed'])==int(x['negative_variance'])+int(x['no_reference_snps'])+int(x['other'] or 0)
 assert abs(float(x['failure_fraction'])-int(x['flagged'])/int(x['receipts']))<1e-12
ins=next(x for x in ledger if x['sleep_trait']=='insomnia')
sha=lambda path:hashlib.sha256(path.read_bytes()).hexdigest()
r={'reference_url':'https://bitbucket.org/nygcresearch/ldetect-data/raw/ac125e47bf7f/EUR/fourier_ls-all.bed','partition_sha256':sha(p/'eur_ld_blocks.bed'),'n_blocks':len(blocks),'per_chromosome':dict(Counter(c for c,s,e in blocks)),'song_exact_block':True,'song_p':nominal,'illustrative_bonferroni_1703':nominal*1703,'threshold_1703':.05/1703,'family_caveat':'Actual Song testing family, overlap correction and build not established by coordinate match. This is an illustration, not a recovered author-adjusted value.','lava_ledger_sha256':sha(p/'lava_gate_ledger.tsv'),'lava_ledger_counts':{'families':len(ledger),'all_gates_failed':True,'receipts':sum(int(x['receipts']) for x in ledger),'negative_variance':sum(int(x['negative_variance']) for x in ledger),'no_reference':sum(int(x['no_reference_snps']) for x in ledger),'other':sum(int(x['other'] or 0) for x in ledger),'insomnia':ins},'lava_caveat':'Saved aggregation only; no fresh per-receipt or matching-locus audit. Flags can overlap; lack of test is not a null result.'}
(p/'audit_receipt.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
