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
import urllib.parse,subprocess,concurrent.futures
comp=str.maketrans('ACGT','TGCA'); items=[]
for rs,(p,ref,alt) in effects.items():
 if len(loc[rs])!=1:continue
 c,pos=next(iter(loc[rs]))
 if c not in ('1','2') or len(seq[rs]['A1'])!=1 or len(seq[rs]['A2'])!=1:continue
 a=next(iter(seq[rs]['A1']));b=next(iter(seq[rs]['A2']))
 if len(a)!=150 or len(b)!=150 or sum(x!=y for x,y in zip(a,b))!=1:continue
 items.append((rs,c,pos,ref,alt,p,a,b))
def one(it):
 rs,c,pos,ref,alt,p,a,b=it
 url='https://api.genome.ucsc.edu/getData/sequence?'+urllib.parse.urlencode({'genome':'hg19','chrom':'chr'+c,'start':pos-80,'end':pos+80}).replace('&',';')
 try:
  d=json.loads(subprocess.check_output(['curl','-fsSL','--retry','2',url],text=True,timeout=25)); dna=d['dna'].upper();hits=[]
  for orientation,sa,sb in [('forward',a,b),('reverse_complement',a.translate(comp)[::-1],b.translate(comp)[::-1])]:
   jvar=next(j for j in range(150) if sa[j]!=sb[j])
   for offset in range(len(dna)-149):
    if pos-1 != (pos-80)+offset+jvar:continue
    reference=dna[offset:offset+150]
    if {sa[jvar],sb[jvar]}!={ref,alt}:continue
    if reference==sa and sb[jvar]!=reference[jvar]:hits.append(orientation)
    if reference==sb and sa[jvar]!=reference[jvar]:hits.append(orientation)
  return ('aligned_'+hits[0] if len(hits)==1 else 'no_or_ambiguous_alignment',int(p<.05))
 except Exception:return ('request_error',None)
with concurrent.futures.ThreadPoolExecutor(max_workers=6) as ex:results=list(ex.map(one,items))
counts=collections.Counter(k for k,v in results);classes=collections.Counter(str(v) for k,v in results if k.startswith('aligned_'))
summary={'examined':len(items),'status_counts':dict(counts),'aligned_adjusted_p_labels':dict(classes),'genome':'UCSC hg19','warning':'Exploratory reference check; does not prove allele-effect direction, independence, phenotype applicability or permission to train'}
open('outputs/stage6_matoba_reference_alignment.json','w').write(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
