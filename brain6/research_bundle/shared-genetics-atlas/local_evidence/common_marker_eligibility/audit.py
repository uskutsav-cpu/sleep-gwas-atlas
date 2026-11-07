import csv,gzip,hashlib,json
from pathlib import Path
p=Path(__file__).resolve().parent; atlas=p.parents[1];base=atlas/'multiscale_embedding/hfrs_extension/heldout_sleep/group_holdout/bounded_model/mvp_feature_transfer/common_support_run';src=base/'MVP.sumstats.gz';prep=base/'prepared_receipt.json';replay=atlas/'local_evidence/song_lava_receipt/cross_trait/marker_overlap/replay/insomnia.json';ver=atlas/'local_evidence/song_lava_receipt/cross_trait/marker_overlap/replay/verification.json'
def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
assert sha(src)==json.loads(prep.read_text())['prepared_sha256']['MVP'];assert sha(replay)==json.loads(ver.read_text())['files_sha256'][str(replay)];r=json.loads(replay.read_text());ids=set(r['snps']);assert len(ids)==r['n_snps']==2741
n=0;hits=set()
with gzip.open(src,'rt',newline='') as f:
 for row in csv.DictReader(f,delimiter='\t'):
  n+=1
  if row['SNP'] in ids:hits.add(row['SNP'])
assert n==991707
(p/'receipt.json').write_text(json.dumps({'mvp_common_rows':n,'replay_ids':len(ids),'matching_rsids':len(hits),'replay_fraction':len(hits)/len(ids),'input_sha256':{str(x.relative_to(atlas)):sha(x) for x in (src,prep,replay,ver)},'matching':'rsID only; no position or allele confirmation','local_result':'not computed'},indent=2)+'\n');print('replay rsIDs',len(ids),'common-input rsIDs',n,'intersection',len(hits))
