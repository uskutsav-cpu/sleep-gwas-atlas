"""Screen two exploratory marker lists against Ye's published FFS COJO leads.
Requires openpyxl, numpy and pysam; input workbook and indexed 1000G reference
are documented in the accompanying review. This is not an HRS association test.
"""
import csv, hashlib, json
from pathlib import Path
import numpy as np
import openpyxl, pysam
root = Path(__file__).resolve().parents[1]
book = root / 'work/ye_supplement.xlsx'
rows = list(openpyxl.load_workbook(book, read_only=True, data_only=True)[' ST1. COJO independent loci'].values)
leads = [dict(chrom=str(r[1]), rsid=str(r[2]), position=int(r[3])) for r in rows[1:] if r[1] is not None]
assert len(leads) == 37
out = {'workbook_sha256': hashlib.sha256(book.read_bytes()).hexdigest(), 'lead_count': len(leads), 'screens': {}}
for filename, chromcol, poscol in [('unvalidated_15_variant_score_prototype.csv','chr','position_GRCh37_inferred'),('unvalidated_22_shared_direction_score.csv','chr_GRCh37_assumed','position_GRCh37_assumed')]:
 with (root/'outputs'/filename).open() as f: variants = list(csv.DictReader(f))
 results=[]
 for v in variants:
  near = [dict(rsid=r['rsid'], distance_bp=r['position']-int(v[poscol])) for r in leads if r['chrom']==v[chromcol] and abs(r['position']-int(v[poscol]))<=1_000_000]
  results.append(dict(variant=v['variant'], exact_id=[r['rsid'] for r in leads if r['rsid']==v['variant']], exact_position=[r['rsid'] for r in leads if r['chrom']==v[chromcol] and r['position']==int(v[poscol])], leads_within_1mb=near))
 out['screens'][filename] = results
panel={r['sample']:r for r in csv.DictReader((root/'work/1kg_phase3_panel.tsv').open(),delimiter='\t')}
url='https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr7.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'
t=pysam.TabixFile(url,index=str(root/'work/stage6_1000g_chr7_phase3.tbi'))
samples=next(x for x in t.header if x.startswith('#CHROM')).split('\t')[9:]
indices=[i for i,s in enumerate(samples) if s in panel and panel[s]['super_pop']=='EUR']; assert len(indices)==503
vectors=[]; alleles=[]
for pos,rs in [(2074122,'rs28528997'),(2083267,'rs62444907')]:
 hits=[x.split('\t') for x in t.fetch('7',pos-1,pos) if x.split('\t')[2]==rs]; assert len(hits)==1
 row=hits[0]; assert len(row[4].split(','))==1
 alleles.append({'rsid':rs,'position':pos,'ref':row[3],'alt':row[4]})
 dos=[]
 for i in indices:
  gt=row[9+i].split(':')[0].replace('|','/').split('/')
  dos.append(sum(int(a)==1 for a in gt) if len(gt)==2 and all(a!='.' for a in gt) else np.nan)
 vectors.append(np.array(dos))
mask=np.isfinite(vectors[0])&np.isfinite(vectors[1]); r=float(np.corrcoef(vectors[0][mask],vectors[1][mask])[0,1])
out['ld']={'source':url,'population':'1000 Genomes phase 3 EUR','n_complete':int(mask.sum()),'sites':alleles,'alt_dosage_r':r,'r_squared':r*r}
(root/'outputs/stage6_ye_ffs_lead_proxy_screen.json').write_text(json.dumps(out,indent=2)+'\n')
print('leads',len(leads),'exact ID matches', [sum(bool(r['exact_id']) for r in x) for x in out['screens'].values()], 'exact position matches', [sum(bool(r['exact_position']) for r in x) for x in out['screens'].values()]); print('LD r²',r*r,'n',int(mask.sum()))
