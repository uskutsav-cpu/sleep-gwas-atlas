import requests,gzip,io,hashlib,csv,json
from pathlib import Path
sources={'insomnia':('https://ndownloader.figshare.com/files/41506827','b7158a866f186395e5b82f4da7208f1d'),'fatigue':('https://ndownloader.figshare.com/files/41506956','e14c0d802d70c5b15458827b63e81315')}
g={}
for label,file in [('score15','outputs/unvalidated_15_variant_score_prototype.csv'),('score22','outputs/unvalidated_22_shared_direction_score.csv')]:
 for r in csv.DictReader(open(file)):
  chrom=r.get('chr',r.get('chr_GRCh37_assumed'));pos=int(r.get('position_GRCh37_inferred',r.get('position_GRCh37_assumed')));pair=frozenset([r['counted_allele'],r['other_allele']]); key=(chrom,pos,pair);g.setdefault(key,{'id':r['variant'],'lists':[],'hits':{}})['lists'].append(label)
class HashReader(io.RawIOBase):
 def __init__(self,raw):self.raw=raw;self.md5=hashlib.md5();self.n=0
 def readable(self):return True
 def readinto(self,b):
  d=self.raw.read(len(b));b[:len(d)]=d;self.md5.update(d);self.n+=len(d);return len(d)
for label,(url,expected) in sources.items():
 r=requests.get(url,stream=True,timeout=(25,180));r.raise_for_status();h=HashReader(r.raw);gz=gzip.GzipFile(fileobj=io.BufferedReader(h));t=io.TextIOWrapper(gz);header=None;pre=[];count=0
 for line in t:
  if line.startswith('##'):
   if line.startswith(('##nCase=','##nControl=','##nCohort=','##cohortList=','##genomeReference=')):pre.append(line.strip())
   continue
  if header is None:header=line.strip().lstrip('#').split('\t');continue
  v=line.rstrip('\n').split('\t');count+=1
  if len(v)!=len(header):raise ValueError((label,count,len(v),len(header)))
  row=dict(zip(header,v));
  if not row['POS'].isdigit():continue
  key=(row.get('CHROM',row.get('CHR')),int(row['POS']),frozenset((row['A1'],row['A2'])))
  if key in g:g[key]['hits'].setdefault(label,[]).append({'a1':row['A1'],'a2':row['A2'],'beta':row['BETA'],'se':row['SE'],'p':row['PVAL'],'fcas':row.get('FCAS'),'fcon':row.get('FCON'),'info':row.get('IMPINFO'),'neff':row.get('NEFF'),'ncase':row.get('NCAS'),'ncontrol':row.get('NCON')})
 assert h.md5.hexdigest()==expected,(label,h.md5.hexdigest());print(label,'records',count,'compressed_bytes',h.n,'md5',h.md5.hexdigest(),'matched',sum(bool(v['hits'].get(label)) for v in g.values()),'headers',pre,flush=True)
 Path('work/mdd_symptoms/'+label+'_metadata.json').write_text(json.dumps({'md5':h.md5.hexdigest(),'rows':count,'header':header,'preheader':pre},indent=2)+'\n')
Path('work/mdd_symptoms/candidate_hits_private.json').write_text(json.dumps([{'chr':k[0],'pos':k[1],'pair':sorted(k[2]),**v} for k,v in g.items()],indent=2)+'\n')
