import csv,gzip,re,collections,openpyxl,json,hashlib
book='work/matoba2026/table3.xlsx';w=openpyxl.load_workbook(book,read_only=True,data_only=True)
s=w['a. Regulatory activity (Veh)'];it=iter(s.iter_rows(values_only=True));h=next(it);loc=collections.defaultdict(set)
for r in it:
 if len(r)>9 and r[1] and r[8] and r[9]:loc[str(r[1])].add((str(r[8]),int(r[9])))
s=w['b. Allelic effect (Veh)'];it=iter(s.iter_rows(values_only=True));h=next(it); effects={}
for r in it:
 if len(r)>7 and r[0] is not None and isinstance(r[7],(float,int)):effects[str(r[0])]=(float(r[7]),str(r[1]),str(r[2]))
seq=collections.defaultdict(lambda:collections.defaultdict(set))
with gzip.open('work/matoba2026/barcode_statistics.gz','rt') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  m=re.search(r'_SNP_(rs\d+)_.*_(A[12])_',r['name'])
  if m and set(r['variant'])<=set('ACGT'):seq[m[1]][m[2]].add(r['variant'])
import urllib.parse,subprocess
comp=str.maketrans('ACGT','TGCA'); answers=[]
for rs,(p,ref,alt) in effects.items():
 if len(loc[rs])!=1:continue
 c,pos=next(iter(loc[rs]));
 if c not in ('1','2') or len(seq[rs]['A1'])!=1 or len(seq[rs]['A2'])!=1:continue
 a=next(iter(seq[rs]['A1']));b=next(iter(seq[rs]['A2']))
 if len(a)!=150 or len(b)!=150 or sum(x!=y for x,y in zip(a,b))!=1:continue
 rec={'snp':rs,'results':{}}
 for genome in ('hg19','hg38'):
  url='https://api.genome.ucsc.edu/getData/sequence?'+urllib.parse.urlencode({'genome':genome,'chrom':'chr'+c,'start':pos-80,'end':pos+80}).replace('&',';')
  data=json.loads(subprocess.check_output(['curl','-fsSL',url],text=True));dna=data['dna'].upper(); hits=[]
  for ori,s in [('forward',a),('reverse_complement',a.translate(comp)[::-1])]:
   for j in range(len(dna)-len(s)+1):
    d=sum(x!=y for x,y in zip(dna[j:j+150],s))
    if d<=1:hits.append((ori,d))
  rec['results'][genome]=hits
 answers.append(rec)
 if len(answers)==6:break
summary={'n':len(answers),'hg19_reverse_complement_with_at_most_one_difference':sum(any(x[0]=='reverse_complement' for x in a['results']['hg19']) for a in answers),'hg38_any_match_with_at_most_one_difference':sum(bool(a['results']['hg38']) for a in answers),'scope':'first six qualifying variants; exploratory; no allele or build guarantee for remaining variants'}
open('outputs/stage6_matoba_orientation_pilot.json','w').write(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
