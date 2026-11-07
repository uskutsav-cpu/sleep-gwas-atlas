import pathlib,csv,hashlib,gzip,json
p=pathlib.Path(__file__).resolve().parent
base=pathlib.Path('/Volumes/Extreme SSD/sleep-frailty-full-autonomous-v5-20261009/frailty_paper/full_autonomous_research_v5/genetics/new_public_validation/raw')
inputs={'insomnia':(base/'finngen_R13_F5_INSOMNIA.gz','sha256','353129dc8252461a2a97087ce6de39548d0b0bf1e488a6ef36da75e9f4690049'),'hfrs':(base/'GCST90707291.tsv','md5','38d9434f50a4493def9d2c9d1cc1b4b3')}
def digest(path,alg):
 h=hashlib.new(alg)
 with path.open('rb') as f:
  while b:=f.read(1048576):h.update(b)
 return h.hexdigest()
report={}
for name,(path,alg,expected) in inputs.items():
 actual=digest(path,alg);assert actual==expected;op=gzip.open if path.suffix=='.gz' else open;dest=p/(name+'_region.tsv');n=0;keep=0
 with op(path,'rt') as f,dest.open('w',newline='') as out:
  head=f.readline().rstrip('\n').split('\t');idx={x:i for i,x in enumerate(head)};chrom='#chrom' if name=='insomnia' else 'chromosome';pos='pos' if name=='insomnia' else 'base_pair_location';w=csv.writer(out,delimiter='\t');w.writerow(head)
  for line in f:
   n+=1;v=line.rstrip('\n').split('\t')
   if v[idx[chrom]]=='7' and 1534487<=int(v[idx[pos]])<=2534487:w.writerow(v);keep+=1
 report[name]={'source_path':str(path),'source_digest':actual,'algorithm':alg,'source_rows':n,'region_rows':keep,'extraction_sha256':digest(dest,'sha256')};print(name,report[name],flush=True)
(p/'gwas_extraction_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
