from pathlib import Path
import csv, json
import numpy as np
import pandas as pd
import pysam
from pyliftover import LiftOver

root = Path(__file__).resolve().parents[1]
d = pd.read_csv(root / 'work/meis1_gtex/meis1_hfrs_aligned.csv')
names = ['rs111736174', 'rs6738144', 'rs12622537', 'rs6546231', 'rs10183293', 'rs62143990']
ix = [d.index[d.rsid == n].tolist() for n in names]
assert all(len(x) == 1 for x in ix)
ix = [x[0] for x in ix]
size = len(d)
mat = {p: np.fromfile(root / f'work/meis1_gtex/meis1_hfrs_{p}_ld.bin', dtype='<f8').reshape(size, size) for p in ['nfe', 'fin']}
with open(root / 'work/1kg_phase3_panel.tsv') as f:
    panel = list(csv.DictReader(f, delimiter='\t'))
groups = {'EUR': {r['sample'] for r in panel if r['super_pop'] == 'EUR'}, 'FIN': {r['sample'] for r in panel if r['pop'] == 'FIN'}}
url = 'https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'
vcf = pysam.VariantFile(url)
lo = LiftOver(str(root / 'work/hg38ToHg19.over.chain.gz'))
genotypes, details = {}, []
for name, i in zip(names, ix):
    row = d.iloc[i]
    conv = lo.convert_coordinate('chr2', int(row.pos38)-1)
    if len(conv) != 1 or conv[0][0] != 'chr2' or conv[0][2] != '+':
        details.append({'rsid':name,'status':'liftover_not_unique_forward'}); continue
    p37 = conv[0][1]+1
    hits = [v for v in vcf.fetch('2',p37-1,p37) if v.pos == p37 and v.id and name in v.id.split(';') and len(v.alleles)==2 and set(v.alleles)=={row.ref,row.alt}]
    if len(hits) != 1:
        details.append({'rsid':name,'status':'no_unique_exact_biallelic_match','hits':len(hits)}); continue
    v=hits[0]; alt_index=v.alleles.index(row.alt); vals={}
    for s in set(vcf.header.samples) & groups['EUR']:
        gt=v.samples[s].get('GT')
        if gt is not None and len(gt)==2 and None not in gt: vals[s]=sum(a==alt_index for a in gt)
    genotypes[name]=vals
    details.append({'rsid':name,'status':'matched','hg19_pos':p37,'vcf_alleles':list(v.alleles),'counted_allele':row.alt,'eur_called':len(vals)})
out=[]
for a in range(len(names)):
    for b in range(a+1,len(names)):
        row={'variant_a':names[a],'variant_b':names[b],'gnomad_nfe_r':float(mat['nfe'][ix[a],ix[b]]),'gnomad_fin_r':float(mat['fin'][ix[a],ix[b]])}
        for group, ids in groups.items():
            common=sorted(ids & set(genotypes.get(names[a],{})) & set(genotypes.get(names[b],{})))
            x=np.array([genotypes[names[a]][s] for s in common]);y=np.array([genotypes[names[b]][s] for s in common])
            row[group+'_n']=len(common)
            row[group+'_r']=float(np.corrcoef(x,y)[0,1]) if len(common)>2 and x.std()>0 and y.std()>0 else None
        out.append(row)
pd.DataFrame(out).to_csv(root/'outputs/stage6_meis1_hfrs_ld_crosscheck.csv',index=False)
(root/'outputs/stage6_meis1_hfrs_ld_crosscheck_receipt.json').write_text(json.dumps({'source':url,'method':'diploid ALT genotype dosage Pearson correlation','variants':details,'gnomad_reference':'previously decoded gnomAD v2.1.1 FIN and NFE matrices','limitation':'different populations and sample sizes; few selected pairs; not a same-reference decoder validation'},indent=2)+'\n')
print(pd.DataFrame(out).to_string(index=False));print(json.dumps(details))
