import csv,hashlib,json,math,re
from pathlib import Path

p=Path(__file__).resolve().parent
receipt=json.loads((p/'prepare_receipt.json').read_text()); n=receipt['common_valid_nonpalindromic']
def digest(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()
for k,path in receipt['output_paths'].items(): assert digest(path)==receipt['output_sha256'][k]
vals={};results=[];deletes={};rounded={}
for physical in ('walking','grip'):
    for sleep in ('questionnaire','registry'):
        name=f'{physical}_{sleep}'; log=p/f'{name}.log'; text=log.read_text()
        assert 'Analysis finished at' in text
        assert int(re.search(r'(\d+) SNPs with valid alleles',text).group(1))==n
        rg,se=map(float,re.search(r'Genetic Correlation: ([\d.-]+) \(([\d.]+)\)',text).groups())
        z,pv=map(float,re.search(r'Z-score: ([\d.-]+)\s+P: ([\d.eE-]+)',text).groups())
        base=p/f'{name}{physical}.sumstats.gz_{sleep}.sumstats.gz'
        v={k:[float(x) for x in Path(str(base)+f'.{k}.delete').read_text().split()] for k in ('hsq1','hsq2','gencov')}
        assert all(len(x)==200 for x in v.values())
        assert all(math.isfinite(x) and x>0 for k in ('hsq1','hsq2') for x in v[k])
        d=[c/math.sqrt(a*b) for a,b,c in zip(v['hsq1'],v['hsq2'],v['gencov'])]
        mean=sum(d)/200; secheck=math.sqrt(199/200*sum((x-mean)**2 for x in d))
        assert abs(secheck-se)<0.00005,(name,secheck,se)
        vals[name]=v; deletes[name]=d;rounded[name]=rg
        results.append({'physical':physical,'sleep':sleep,'rg':rg,'se':se,'logged_p':pv,
                        'bonferroni_four_from_rounded_p':min(1,pv*4),'retained':n,'log_sha256':digest(log)})
for sleep in ('questionnaire','registry'):
    assert vals[f'walking_{sleep}']['hsq2']==vals[f'grip_{sleep}']['hsq2']
contrasts=[]
for physical in ('walking','grip'):
    q=f'{physical}_questionnaire'; r=f'{physical}_registry'
    assert vals[q]['hsq1']==vals[r]['hsq1']
    d=[a-b for a,b in zip(deletes[q],deletes[r])];avg=sum(d)/200
    se=math.sqrt(199/200*sum((x-avg)**2 for x in d)); delta=rounded[q]-rounded[r]
    pv=math.erfc(abs(delta/se)/math.sqrt(2))
    contrasts.append({'physical':physical,'difference_questionnaire_minus_registry':delta,
                      'paired_jackknife_se':se,'ci95_unadjusted':[delta-1.95996398454*se,delta+1.95996398454*se],
                      'p':pv,'bonferroni_two':min(1,2*pv),'identical_physical_heritability_deletes':True})
with (p/'rg_deletes.tsv').open('w',newline='') as f:
    w=csv.writer(f,delimiter='\t');w.writerow(['block']+list(deletes))
    for i in range(200):w.writerow([i+1]+[deletes[k][i] for k in deletes])
out={'correlations':results,'contrasts':contrasts,'source_prepare_receipt_sha256':digest(p/'prepare_receipt.json'),
     'common_markers':n,'deletion_blocks':200,'note':'Point rg and p are read from rounded log output; selection, phenotype and source differences limit inference.'}
(p/'results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
