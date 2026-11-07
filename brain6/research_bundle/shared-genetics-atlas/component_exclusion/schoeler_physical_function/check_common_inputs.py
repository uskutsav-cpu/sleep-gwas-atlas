import csv, gzip, hashlib, json
from pathlib import Path
here = Path(__file__).resolve().parent
atlas = here.parents[1]
base = atlas/'multiscale_embedding/hfrs_extension/heldout_sleep/group_holdout/bounded_model/mvp_feature_transfer/common_support_run'
ref = {}
result = {}
complement = str.maketrans('ACGT','TGCA')
for name in ('MVP','FI','FFS','HFRS'):
    path = base/(name+'.sumstats.gz')
    sha=hashlib.sha256(path.read_bytes()).hexdigest()
    with gzip.open(path,'rt',newline='') as f:
        reader=csv.DictReader(f,delimiter='\t')
        seen={}
        for row in reader:
            snp=row['SNP']; assert snp not in seen
            a,b=row['A1'],row['A2']; assert len(a)==len(b)==1 and a in 'ACGT' and b in 'ACGT' and a!=b
            pair=''.join(sorted((a,b)))
            assert pair not in ('AT','CG')
            seen[snp]=pair
    if not ref: ref=seen
    missing=set(ref)-set(seen);extra=set(seen)-set(ref)
    incompatible=sum(not (pair==ref[snp] or ''.join(sorted(pair.translate(complement)))==ref[snp]) for snp,pair in seen.items() if snp in ref)
    assert not missing and not extra and not incompatible and len(seen)==991707
    result[name]={'sha256':sha,'unique_snp_count':len(seen),'missing_vs_MVP':len(missing),'extra_vs_MVP':len(extra),'allele_pair_incompatible':incompatible}
(here/'common_inputs_receipt.json').write_text(json.dumps({'inputs':result,'interpretation':'The prepared four inputs have identical SNP sets and compatible non-palindromic allele pairs. The streamed source SNP list was not retained; LDSC logs additionally report equal allele-valid counts, so the exact output marker sets were not independently exported.'},indent=2)+'\n')
print(json.dumps(result,indent=2))
