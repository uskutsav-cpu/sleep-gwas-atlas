import csv,json,hashlib,pathlib
import numpy as np,pysam
rows=[x for x in csv.DictReader(open('outputs/stage5_smad5_splice_snp_ld_reference_coverage.csv')) if x['status']=='unique allele match'];selected=[rows[i*(len(rows)-1)//9] for i in range(10)];fam=[x.split()[1] for x in open('work/stage4_magma_reference/g1000_eur.fam')];n=len(fam);step=(n+3)//4;bed=np.memmap('work/stage4_magma_reference/g1000_eur.bed',mode='r',dtype='uint8');lookup=np.array([2,np.nan,1,0]);t=pysam.TabixFile('https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr5.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz',index='work/stage5_1000g_chr5_phase3.tbi');header=next(x for x in t.header if x.startswith('#CHROM')).split('\t')[9:];indices=[header.index(x) for x in fam];results=[]
for r in selected:
 pos=int(r['pos37']);expected={r['bim_a1'],r['bim_a2']};vs=[]
 for line in t.fetch('5',pos-1,pos):
  z=line.split('\t');alts=z[4].split(',')
  if int(z[1])==pos and len(alts)==1 and {z[3],alts[0]}==expected:vs.append(z)
 if len(vs)!=1:results.append({'variant_GRCh38':r['variant_GRCh38'],'status':'no unique biallelic match','compared':0,'discordant':''});continue
 z=vs[0];target=r['alt38'];allele=0 if target==z[3] else 1;gtcol=z[8].split(':').index('GT');dos=[]
 for i in indices:
  gt=z[9+i].split(':')[gtcol].replace('|','/').split('/');dos.append(sum(int(v)==allele for v in gt) if len(gt)==2 and all(v.isdigit() for v in gt) else np.nan)
 offset=3+int(r['bim_index_zero_based'])*step;raw=bed[offset:offset+step];code=((raw[:,None]>>np.array([0,2,4,6],dtype=np.uint8))&3).reshape(-1)[:n];a=lookup[code];a=2-a if r['bim_a1']!=target else a;b=np.array(dos);ok=np.isfinite(a)&np.isfinite(b);results.append({'variant_GRCh38':r['variant_GRCh38'],'status':'compared','compared':int(ok.sum()),'discordant':int(np.sum(a[ok]!=b[ok]))})
with open('outputs/stage5_smad5_ld_reference_genotype_crosscheck.csv','w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(results[0]));w.writeheader();w.writerows(results)
print(json.dumps({'snvs_sampled':len(selected),'snvs_compared':sum(x['status']=='compared' for x in results),'genotype_comparisons':sum(x['compared'] for x in results),'discordant':sum(x['discordant'] for x in results if isinstance(x['discordant'],int))},indent=2))
