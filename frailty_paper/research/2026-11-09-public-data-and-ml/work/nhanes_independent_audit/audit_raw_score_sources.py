"""Compare supplied score alleles and weights with original GWAS inputs.

Checks exact stated rsID, chromosome, coordinate and unordered allele pair.
FFS A1 interpretation requires the separate published COJO audit.
"""
import csv, gzip, hashlib, json, math
from pathlib import Path

root=Path(__file__).resolve().parent
package=root/'finalized_study'
raw=Path('/Users/swethasunilkumar/Documents/Codex/2026-09-19/can/work/sleep-gwas-atlas-v03/data/raw')
ffs=Path('/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/frailty_paper/data/gwas/fried_frailty_score/Frailty_GWAS_maf0.01_info0.9_hwe1e-6_FUMA.tsv')
scores=[]
for n in (15,22):
    with (package/f'score_{n}_variants_standardized.csv').open() as f:
        scores.extend(dict(row,model=n) for row in csv.DictReader(f))
ids={r['rsid'] for r in scores}
sources={
 'FI':(raw/'frailty.txt.gz', {'id':'variant_id','chr':'chromosome','pos':'base_pair_location','a1':'effect_allele','a2':'other_allele','weight':'beta'}),
 'insomnia':(raw/'insomnia.txt.gz', {'id':'SNP','chr':'CHR','pos':'BP','a1':'A1','a2':'A2','weight':'OR'}),
 'FFS':(ffs, {'id':'SNP','chr':'CHROM','pos':'POS','a1':'A1','a2':'A2','weight':'BETA'})}
weights={'FI':'FI_beta_weight','insomnia':'insomnia_beta_weight','FFS':'FFS_beta_weight_inferred'}
out=[];receipt=[]
for trait,(path,cols) in sources.items():
    digest=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(4*1024*1024),b''):digest.update(block)
    opening=gzip.open if path.suffix=='.gz' else open
    hits={};n=0
    with opening(path,'rt',encoding='utf8',newline='') as f:
        reader=csv.DictReader(f,delimiter='\t')
        assert set(cols.values())<=set(reader.fieldnames),(path,reader.fieldnames)
        for r in reader:
            n+=1
            if r[cols['id']] in ids:hits.setdefault(r[cols['id']],[]).append(r)
    receipt.append({'trait':trait,'path':str(path),'sha256':digest.hexdigest(),'rows':n})
    for row in scores:
        possible=[r for r in hits.get(row['rsid'],[]) if r[cols['chr']].removeprefix('chr')==row['chromosome'] and r[cols['pos']]==row['position']]
        result={'model':row['model'],'rsid':row['rsid'],'trait':trait,'coordinate_matches':len(possible),'allele_alignment':'','source_value':'','source_effect':'','score_weight':row[weights[trait]],'difference':'','status':'not_unique'}
        if len(possible)==1:
            r=possible[0];a=r[cols['a1']].upper();b=r[cols['a2']].upper();c=row['effect_counted_allele'].upper();d=row['other_allele'].upper()
            orient=1 if (a,b)==(c,d) else -1 if (a,b)==(d,c) else 0
            value=float(r[cols['weight']]); beta=math.log(value) if trait=='insomnia' else value
            diff=abs(orient*beta-float(row[weights[trait]])) if orient else None
            result.update(allele_alignment='same' if orient==1 else 'swapped' if orient==-1 else 'mismatch',source_value=value,source_effect=a,difference=diff if diff is not None else '',status='match' if diff is not None and diff<1e-10 else 'difference_or_mismatch')
        out.append(result)
with (root/'raw_score_source_audit.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=out[0]);w.writeheader();w.writerows(out)
(root/'raw_score_source_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
from collections import Counter
print(Counter((r['trait'],r['status']) for r in out))
