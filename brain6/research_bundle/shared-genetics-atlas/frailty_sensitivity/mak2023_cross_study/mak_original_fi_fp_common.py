import csv,gzip,math,pathlib,collections,json,hashlib
mak=pathlib.Path('/Volumes/Extreme SSD/Codex Archive/mak2023_same_cohort');old=pathlib.Path('/Volumes/Extreme SSD/Codex Archive/hfrs_fi_ffs_common_support');out=mak/'fourway';out.mkdir(exist_ok=True)
old1=gzip.open(old/'FI.sumstats.gz','rt');old2=gzip.open(old/'FFS.sumstats.gz','rt');m1=gzip.open(mak/'mak_fi_common.sumstats.gz','rt');m2=gzip.open(mak/'mak_fp_common.sumstats.gz','rt'); oldrows={};dups=set()
for a,b in zip(csv.DictReader(old1,delimiter='\t'),csv.DictReader(old2,delimiter='\t')):
 assert a['SNP']==b['SNP'] and (a['A1'],a['A2'])==(b['A1'],b['A2'])
 assert a['SNP'] not in oldrows
 oldrows[a['SNP']]=(a,b)
outputs=[gzip.open(out/(n+'.sumstats.gz'),'wt') for n in ('old_FI','old_FFS','mak_FI','mak_FP')]
for f in outputs:f.write('SNP\tA1\tA2\tN\tZ\n')
n=0;bad=0;complement=str.maketrans('ACGT','TGCA')
for x,y in zip(csv.DictReader(m1,delimiter='\t'),csv.DictReader(m2,delimiter='\t')):
 assert x['SNP']==y['SNP'] and (x['A1'],x['A2'])==(y['A1'],y['A2']);oldpair=oldrows.get(x['SNP'])
 if oldpair is None:continue
 a,b=oldpair;aa,bb=x['A1'],x['A2'];oa,ob=a['A1'],a['A2']
 if (oa,ob)==(aa,bb) or (oa.translate(complement),ob.translate(complement))==(aa,bb):sign=1
 elif (ob,oa)==(aa,bb) or (ob.translate(complement),oa.translate(complement))==(aa,bb):sign=-1
 else:bad+=1;continue
 for f,z in zip(outputs,(a,b,x,y)):
  zz=float(z['Z'])*(sign if z is a or z is b else 1);f.write(f"{x['SNP']}\t{aa}\t{bb}\t{z['N']}\t{zz:.10g}\n")
 n+=1
for f in outputs+[old1,old2,m1,m2]:f.close()
assert n>500000
meta={'n_common':n,'allele_mismatches_excluded':bad,'description':'Four-source common SNPs previously restricted by HFRS and FinnGen insomnia availability. A1/A2 aligned to Mak, including strand complements at nonpalindromic positions; old source z scores already rounded.'}
(out/'support.json').write_text(json.dumps(meta,indent=2)+'\n');print(json.dumps(meta,indent=2))
