import csv,json,collections,math
from pathlib import Path
import numpy as np,pysam
from pyliftover import LiftOver
base=Path('work/meis1_gtex'); sets={}
for tag in ('FI','FFS','Jansen'):
 with (base/f'{tag}_eqtl_wide_matched.tsv').open() as f:
  rows=list(csv.DictReader(f,delimiter='\t'))
 sets[tag]={(int(r['pos38']),r['ref38'],r['alt38']):r for r in rows if 66400000<=int(r['pos38'])<=66580000 and {r['ref38'],r['alt38']} not in ({'A','T'},{'C','G'})}
common=sorted(set.intersection(*(set(v) for v in sets.values())))
q={}
for line in (base/'QTD000166_MEIS1AS3_65400000_67500000.tsv').read_text().splitlines():
 x=line.split('\t');q[(int(x[2]),x[3],x[4])]=x
lo=LiftOver('work/hg38ToHg19.over.chain.gz');u='https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz';f=pysam.VariantFile(u)
with open('work/1kg_phase3_panel.tsv') as h: ids={r['sample'] for r in csv.DictReader(h,delimiter='\t') if r['super_pop']=='EUR'}
ids=sorted(ids & set(f.header.samples));by37={}
for k in common:
 m=lo.convert_coordinate('chr2',k[0]-1)
 if len(m)==1 and m[0][0]=='chr2' and m[0][2]=='+':by37[(m[0][1]+1,k[1],k[2])]=k
xs={}
for x in f.fetch('2',66627129,66807134):
 if len(x.alleles)!=2:continue
 k=by37.get((x.pos,x.alleles[0],x.alleles[1]))
 if not k:continue
 g=[]
 for s in ids:
  gt=x.samples[s].get('GT');g.append(sum(z==1 for z in gt) if gt and len(gt)==2 and None not in gt else float('nan'))
 if np.all(np.isfinite(g)) and np.std(g)>0:xs[k]=np.array(g,dtype=float)
out=[];mat=[];excluded=collections.Counter()
for k in common:
 if k not in xs:excluded['reference_missing_or_mismatch']+=1;continue
 x=q[k];af=float(x[12])/float(x[13]);afref=float(xs[k].mean()/2)
 if abs(af-afref)>0.2:excluded['qtl_reference_af_difference_gt_0.2']+=1;continue
 row={'id':f'2:{k[0]}:{k[1]}:{k[2]}','rsid':x[-1],'pos38':k[0],'ref':k[1],'alt':k[2],'qtl_alt_af':af,'ref_alt_af':afref,'qtl_beta':float(x[9]),'qtl_se':float(x[10])}
 try:
  for tag in sets:
   r=sets[tag][k];e=float(r['gwas_effect']);se=float(r['gwas_se']);beta=math.log(e) if tag=='Jansen' else e
   if r['gwas_effect_allele']!=k[2]:beta=-beta
   if not math.isfinite(beta) or not math.isfinite(se) or se<=0:raise ValueError('bad stats')
   row[tag+'_beta']=beta;row[tag+'_se']=se
 except (ValueError,KeyError):excluded['invalid_statistics']+=1;continue
 out.append(row);mat.append(xs[k])
arr=np.corrcoef(np.stack(mat));np.asarray(arr,dtype='<f8').tofile(base/'meis1_region_eur_ld.bin')
with (base/'meis1_susie_region_aligned.csv').open('w',newline='') as h:w=csv.DictWriter(h,fieldnames=out[0].keys());w.writeheader();w.writerows(out)
res={'region_GRCh38':'chr2:66400000-66580000','common_nonpalindromic':len(common),'reference_matched_after_qc':len(out),'n_reference':len(ids),'excluded':dict(excluded),'qtl_reference_max_af_difference':max(abs(r['qtl_alt_af']-r['ref_alt_af']) for r in out),'reference_ld_rank':int(np.linalg.matrix_rank(arr)),'eigen_min':float(np.linalg.eigvalsh(arr).min()),'caveat':'postselected narrow region; 1000G EUR external LD, qtl vs ref AF only; GWAS AF and allele source representations unverified; not score validation'}
Path('outputs/stage6_meis1_susie_region_preparation.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res,indent=2))
