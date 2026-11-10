#!/usr/bin/env python3
"""Extract exact GRCh37 candidate records only after full checksum verification."""
import csv,gzip,hashlib,json
from pathlib import Path

root=Path(__file__).resolve().parent.parent; out=root/'outputs';source=root/'work/GCST90319502.partial.snapshot.gz'
expected='0fc2688cf786984fb5f147f6eb52628a'
h=hashlib.md5()
with source.open('rb') as f:
    for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
if h.hexdigest()!=expected:
    raise SystemExit(f'Source incomplete or checksum mismatch: {h.hexdigest()}; no outputs created')
def rows(path):
    with path.open(newline='') as f:return list(csv.DictReader(f))
specs=[('prototype15','unvalidated_15_variant_score_prototype.csv','chr','position_GRCh37_inferred'),('shared22','unvalidated_22_shared_direction_score.csv','chr_GRCh37_assumed','position_GRCh37_assumed')]
targets=[];by_position={}
for label,file,ch,pos in specs:
    for r in rows(out/file):
        target={'score':label,'variant':r['variant'],'chromosome':r[ch],'position':int(r[pos]),'counted_allele':r['counted_allele'],'other_allele':r['other_allele'],'position_rows':0,'exact_rows':0,'exact_records':[]}
        targets.append(target);by_position.setdefault((target['chromosome'],target['position']),[]).append(target)
assert len(targets)==37
# The paper's published lead provides a source-specific orientation control.
control={'score':'paper_control','variant':'rs10749438','chromosome':'10','position':123810832,'counted_allele':'A','other_allele':'G','position_rows':0,'exact_rows':0,'exact_records':[]}
targets.append(control);by_position.setdefault(('10',123810832),[]).append(control)
required={'chromosome','base_pair_location','effect_allele','other_allele','beta','standard_error','effect_allele_frequency','p_value','variant_id','n'}
scanned=0
with gzip.open(source,'rt',newline='') as f:
    reader=csv.DictReader(f,delimiter='\t');assert required.issubset(reader.fieldnames),reader.fieldnames
    for row in reader:
        scanned+=1
        try:key=(row['chromosome'].removeprefix('chr'),int(row['base_pair_location']))
        except (ValueError,KeyError):continue
        for t in by_position.get(key,[]):
            t['position_rows']+=1
            if {row['effect_allele'].upper(),row['other_allele'].upper()}=={t['counted_allele'],t['other_allele']}:
                t['exact_rows']+=1;t['exact_records'].append({k:row[k] for k in required})
output=[]
for t in targets:
    row={k:v for k,v in t.items() if k!='exact_records'}
    if t['exact_rows']==1:
        v=t['exact_records'][0];row.update({f'source_{k}':v[k] for k in required});row['beta_counted_if_source_coding_valid']=float(v['beta'])*(1 if v['effect_allele']==t['counted_allele'] else -1)
    output.append(row)
fields=['score','variant','chromosome','position','counted_allele','other_allele','position_rows','exact_rows',*[f'source_{k}' for k in sorted(required)],'beta_counted_if_source_coding_valid']
with (out/'stage6_japanese_knee_candidate_matches.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(output)
ctrl=next(r for r in output if r['score']=='paper_control')
result={'source':'https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90319001-GCST90320000/GCST90319502/GCST90319502.tsv.gz','observed_md5':h.hexdigest(),'metadata_expected_md5':expected,'checksum_matches':True,'rows_scanned':scanned,'exact_single_match_by_score':{label:sum(r['score']==label and r['exact_rows']==1 for r in output) for label,*_ in specs},'paper_control':ctrl,'limitations':'Unharmonized source; exact position/allele set alone does not establish strand, multiallelic contrast, independence, or score performance. Beta orientation assumes source effect-allele convention. Compare control with the published ALT-A beta −0.148, SE 0.027 and p 4.24e−8 before interpreting target effects.'}
(out/'stage6_japanese_knee_candidate_receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
