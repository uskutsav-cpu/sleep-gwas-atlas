import csv,hashlib,json,pathlib
p=pathlib.Path('outputs/shared-genetics-atlas');out=p/'placo_diagnostic/completed_audit_fi_only/prior_literature';x=p/'placo_diagnostic/completed_audit_fi_only/candidate_positions.tsv';y=p/'evidence/song124_mvp_hfrs_joint_screen.csv'
with x.open() as f:a=list(csv.DictReader(f,delimiter='\t'))
with y.open() as f:b=list(csv.DictReader(f))
assert len(a)==16 and len(b)==124 and len({r['variant'] for r in b})==124
rows=[]
for r in a:
 same=[q for q in b if q['variant']==r['SNP']]
 chrom=[q for q in b if q['song_chr']==r['CHR'] and q['song_position_reported'].isdigit()]
 nearest=min(chrom,key=lambda q:abs(int(q['song_position_reported'])-int(r['BP_GRCH37']))) if chrom else None
 rows.append({'SNP':r['SNP'],'chr37':r['CHR'],'bp37':r['BP_GRCH37'],'exact_song_rsid':same[0]['variant'] if len(same)==1 else '', 'exact_song_rank':same[0]['song_cpassoc_rank'] if len(same)==1 else '', 'nearest_song_rsid_same_chr':nearest['variant'] if nearest else '', 'nearest_song_reported_position':nearest['song_position_reported'] if nearest else '', 'absolute_distance_bp_assuming_song_grch37':abs(int(nearest['song_position_reported'])-int(r['BP_GRCH37'])) if nearest else ''})
with (out/'song_crosswalk.tsv').open('w',newline='') as f:w=csv.DictWriter(f,list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
(out/'receipt.json').write_text(json.dumps({'candidate_input_sha256':hashlib.sha256(x.read_bytes()).hexdigest(),'song_lookup_sha256':hashlib.sha256(y.read_bytes()).hexdigest(),'candidate_n':len(a),'song_n':len(b),'exact_rsid_matches':sum(bool(r['exact_song_rsid']) for r in rows),'caveat':'Song build remains assumed, nearest position is not LD nor independent replication; not a comprehensive literature search.'},indent=2)+'\n')
print('exact',[(r['SNP'],r['exact_song_rank']) for r in rows if r['exact_song_rsid']])
