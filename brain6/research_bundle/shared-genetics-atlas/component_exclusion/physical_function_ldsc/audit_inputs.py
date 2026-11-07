import csv, gzip, hashlib, json, math
from pathlib import Path

HERE = Path(__file__).resolve().parent
paths = {'insomnia': Path('/Volumes/Extreme SSD/Codex Archive/hfrs_ukb_sensitivity/S.sumstats.gz'),
         'walking': Path('/Volumes/Brain6Work/neale_924_walking.ldsc.bgz'),
         'grip': Path('/Volumes/Brain6Work/neale_46irnt_grip.ldsc.bgz')}
def digest(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''): h.update(chunk)
    return h.hexdigest()
def read(p):
    seen={}; dup=invalid=0; ns=set(); count=0
    with gzip.open(p,'rt') as f:
        reader=csv.DictReader(f,delimiter='\t')
        assert reader.fieldnames == ['SNP','A1','A2','N','Z'], reader.fieldnames
        for x in reader:
            count+=1; a,b=x['A1'],x['A2']
            if a not in 'ACGT' or b not in 'ACGT' or a==b or len(a)!=1 or len(b)!=1 or not math.isfinite(float(x['Z'])) or not (math.isfinite(float(x['N'])) and float(x['N'])>0): invalid+=1
            if x['SNP'] in seen: dup+=1
            seen[x['SNP']]=(a,b); ns.add(float(x['N']))
    return seen,{'rows':count,'unique_ids':len(seen),'duplicate_ids':dup,'invalid_rows':invalid,'unique_n_count':len(ns),'n_min':min(ns),'n_max':max(ns),'sha256':digest(p),'path':str(p)}
ref,r=read(paths['insomnia']); out={'insomnia':r}
comp=str.maketrans('ACGT','TGCA')
for k in ('walking','grip'):
    d,info=read(paths[k]); shared=set(d)&set(ref); match=ambiguous=0
    for snp in shared:
        a,b=d[snp]; c,e=ref[snp]
        if {a,b} == {c,e} or {a.translate(comp),b.translate(comp)} == {c,e}: match+=1
        if {a,b} in ({'A','T'},{'C','G'}): ambiguous+=1
    info.update(shared_rsid=len(shared),allele_pair_or_strand_matches=match,ambiguous_shared=ambiguous)
    out[k]=info
    assert not info['duplicate_ids'] and not info['invalid_rows'] and match>500000
(HERE/'input_receipt.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
