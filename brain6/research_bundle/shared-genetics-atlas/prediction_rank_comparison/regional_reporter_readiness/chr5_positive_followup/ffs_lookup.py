import csv,hashlib,json,pathlib
p=pathlib.Path(__file__).resolve().parent;source=pathlib.Path('/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/frailty_paper/data/gwas/fried_frailty_score/Frailty_GWAS_maf0.01_info0.9_hwe1e-6_FUMA.tsv');expected='912e290b2064999c7147cbe3597f87cbcb14b25d387a8b6db2f25862e7274e55'
h=hashlib.sha256();hits=[];n=0
with source.open('rb') as f:
 header=f.readline();h.update(header);cols=header.decode().rstrip('\n\r').split('\t');ci=cols.index('CHROM');pi=cols.index('POS');si=cols.index('SNP')
 for line in f:
  h.update(line);n+=1;z=line.decode().rstrip('\n\r').split('\t')
  if z[si]=='rs12659965' or (z[ci]=='5' and z[pi]=='103904399'):hits.append(dict(zip(cols,z)))
assert h.hexdigest()==expected and n==8883488
(p/'ffs_receipt.json').write_text(json.dumps(dict(source=str(source),sha256=h.hexdigest(),rows_scanned=n,matches=hits),indent=2)+'\n');print('source hash verified; rows',n,'matches',len(hits))
