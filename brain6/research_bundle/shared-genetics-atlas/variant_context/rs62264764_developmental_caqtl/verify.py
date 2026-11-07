import csv,hashlib,json
from pathlib import Path
p=Path(__file__).resolve().parent
r=json.loads((p/'scan_receipt.json').read_text());q=json.loads((p/'peak_receipt.json').read_text());src=json.loads((p/'retrieval_receipt.json').read_text())
assert r['total_candidate_rows']==0 and all(v['full_gzip_read'] and v['inner_sha256_verified'] for v in r['datasets'].values())
assert all(v['crc_verified'] for v in src['members']) and not src['archive_fully_hashed']
with (p/'overlap_peaks.tsv').open() as f: rows=list(csv.DictReader(f,delimiter='\t'))
assert len(rows)==2 and {x['cell'] for x in rows}=={'neuron','progenitor'}
for x in rows:
 start,end=int(x['start']),int(x['end']);v=117920728;window=(v-251,v+250)
 assert x['peak']=='chr3_117920181_117920630' and start<window[1] and end>window[0] and int(x['source_rows'])>0
assert hashlib.sha256((p/'candidate_rows.tsv').read_bytes()).hexdigest()==r['candidate_rows_sha256']
assert hashlib.sha256((p/'overlap_peaks.tsv').read_bytes()).hexdigest()==q['file_sha256']
result={'candidate_rows':0,'peak_rows':len(rows),'predicted_window_bed_assumption':[117920477,117920978],'peak_intersection_bases_under_bed_assumption':153,'gzip_complete_and_member_crcs_verified':True,'whole_archive_verified':False,'scope':'These checks do not establish variant effects, no testing, causality, or full source build and peak-coordinate semantics.'}
(p/'verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
