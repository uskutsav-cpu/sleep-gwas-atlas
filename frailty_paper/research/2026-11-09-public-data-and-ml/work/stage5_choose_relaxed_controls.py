import csv,collections,heapq,json,math,numpy as np
markers=[]
for x in csv.DictReader(open('outputs/shared22_eight_source_evidence_matrix.csv')):markers.append({'candidate':x['variant'],'chr':x['chr_GRCh37_assumed'],'pos':int(x['position_GRCh37_assumed']),'alleles':{x['score_counted_allele'],x['other_allele']},'window':[],'anchors':[]})
x=next(x for x in csv.DictReader(open('outputs/integrated_15_candidate_review.csv')) if x['variant']=='rs323509');markers.append({'candidate':'rs323509','chr':x['GRCh37_chr'],'pos':int(x['GRCh37_pos']),'alleles':{x['counted_allele'],x['other_allele']},'window':[],'anchors':[]})
markers=[x for x in markers if x['candidate'] in ('rs12536335','rs6452792')]
by=collections.defaultdict(list)
for x in markers:by[x['chr']].append(x)
base='work/stage4_magma_reference/g1000_eur'
with open(base+'.bim') as f:
 for i,line in enumerate(f):
  c,rs,_,p,a,b=line.split();p=int(p)
  if len(a)!=1 or len(b)!=1 or a not in 'ACGT' or b not in 'ACGT' or not rs.startswith('rs'):continue
  for x in by.get(c,[]):
   if abs(p-x['pos'])<=1000000:
    item=(i,rs,p,a,b);x['window'].append(item)
    if p==x['pos'] and {a,b}==x['alleles']:x['anchors'].append(item)
variants=i+1;n=503;stride=126;bed=np.memmap(base+'.bed',mode='r',dtype='uint8');assert len(bed)==3+variants*stride and bytes(bed[:3])==bytes([0x6c,0x1b,1]);look=np.array([2,np.nan,1,0],dtype='float32')
def decode(ix):
 ix=np.asarray(ix,dtype='int64');b=bed[3+ix[:,None]*stride+np.arange(stride)[None,:]];return look[((b[:,:,None]>>np.array([0,2,4,6],dtype='uint8'))&3).reshape(len(ix),-1)[:,:503]]
pool={};ids=set()
for x in markers:
 assert len(x['anchors'])==1,(x['candidate'],len(x['anchors']));anchor=x['anchors'][0];y=decode([anchor[0]])[0];am=min(float(np.nanmean(y))/2,1-float(np.nanmean(y))/2);options=[]
 for j in range(0,len(x['window']),1024):
  chunk=x['window'][j:j+1024];g=decode([z[0] for z in chunk]);valid=np.isfinite(g)&np.isfinite(y)[None,:];nn=valid.sum(1);gx=np.where(valid,g,0);gy=np.where(valid,y[None,:],0);sx=gx.sum(1);sy=gy.sum(1);cov=(gx*gy).sum(1)-sx*sy/np.maximum(nn,1);vx=(gx*gx).sum(1)-sx*sx/np.maximum(nn,1);vy=(gy*gy).sum(1)-sy*sy/np.maximum(nn,1);den=vx*vy;r2=np.divide(cov*cov,den,out=np.full_like(cov,np.inf),where=den>0);fr=np.divide(sx,2*nn,out=np.zeros_like(sx),where=nn>0);maf=np.minimum(fr,1-fr)
  for k in np.flatnonzero((nn>=400)&(r2<=.2)&(maf>=.01)&(abs(maf-am)<=.05)):
   idx,rs,p,a,b=chunk[int(k)];dist=abs(p-x['pos'])
   if dist<10000 or rs==x['candidate']:continue
   options.append({'candidate':x['candidate'],'chr37':x['chr'],'anchor_pos37':x['pos'],'anchor_maf_reference':am,'control_rsid':rs,'control_pos37':p,'control_bim_alleles':a+'/'+b,'control_maf_reference':float(maf[k]),'r2_to_candidate':float(r2[k]),'n_pairwise':int(nn[k]),'distance_bp':dist,'maf_diff':abs(float(maf[k])-am)})
 options.sort(key=lambda z:(z['maf_diff'],z['distance_bp'],z['control_rsid']));pool[x['candidate']]=options[:500];ids.update(z['control_rsid'] for z in options[:500]);print(x['candidate'],'pool',len(options),'shortlist',len(pool[x['candidate']]),flush=True)
traits={}
for trait in ['atkins_fi','ye_ffs','jansen_insomnia']:
 vals={};dup=set()
 with open('work/stage4_magma_reference/'+trait+'.pval') as f:
  f.readline()
  for line in f:
   x=line.split()
   if x and x[0] in ids:
    if x[0] in vals:dup.add(x[0])
    else:vals[x[0]]=float(x[1])
 for x in dup:vals.pop(x,None)
 traits[trait]=vals;print(trait,'matched',len(vals),'duplicate_ids_excluded',len(dup),flush=True)
selected=[];used=set()
for x in markers:
 opts=[z for z in pool[x['candidate']] if z['control_rsid'] not in used and all(z['control_rsid'] in vals and vals[z['control_rsid']]>.05 for vals in traits.values())]
 if not opts:
  selected.append({'candidate':x['candidate'],'status':'no qualifying control in top 500; unresolved'});continue
 z=opts[0];used.add(z['control_rsid']);selected.append({**z,**{t+'_p':str(vals[z['control_rsid']]) for t,vals in traits.items()},'status':'candidate control; GRCh38 mapping and reference check pending','interpretation':'wider-window, relaxed-LD frequency-matched comparator with p>0.05 in 3 inputs, not a proven neutral variant'})
fields=['candidate','chr37','anchor_pos37','anchor_maf_reference','control_rsid','control_pos37','control_bim_alleles','control_maf_reference','r2_to_candidate','n_pairwise','distance_bp','maf_diff','atkins_fi_p','ye_ffs_p','jansen_insomnia_p','status','interpretation']
with open('outputs/stage5_relaxed_control_candidates_grch37.csv','w',newline='') as f:w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(selected)
print('selected',len(used),'unresolved',len(markers)-len(used),flush=True)
