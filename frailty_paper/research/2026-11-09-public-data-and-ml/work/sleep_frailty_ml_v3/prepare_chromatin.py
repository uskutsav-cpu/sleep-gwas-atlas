"""Annotate public physical reference chromatin, never reporter target readouts.

GSE113480 is hg19; ENCODE SK-N-SH is GRCh38. Coordinate availability and
peak absence are distinct. Profiles are covariates, not new labelled subjects.
"""
from pathlib import Path
import sys,csv,json,gzip,hashlib,collections,bisect,math
import numpy as np
ROOT=Path(__file__).resolve().parents[2];WORK=Path(__file__).resolve().parent;OUT=ROOT/'outputs/sleep_frailty_ai_v3'
sys.path.insert(0,str(ROOT/'work/stage5_env/lib/python3.13/site-packages'));from pyliftover import LiftOver
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load_bed(path,has_signal):
    d=collections.defaultdict(list)
    with gzip.open(path,'rt') as f:
        for line in f:
            if line.startswith(('#','track','browser')):continue
            row=line.rstrip().split('\t');s,e=int(row[1]),int(row[2]);assert 0<=s<e
            value=float(row[6]) if has_signal else np.nan
            assert not has_signal or np.isfinite(value) and value>=0
            d[row[0]].append((s,e,value))
    out={}
    for c,rr in d.items():
        a=np.asarray(sorted(rr),dtype='float64');s=a[:,0].astype('int64');e=a[:,1].astype('int64');out[c]=(s,e,np.maximum.accumulate(e),a[:,2])
    return out
def candidates(index,start,end):
    s,e,prefix,v=index;lo=np.searchsorted(prefix,start,side='right');hi=np.searchsorted(s,end,side='left');ix=np.arange(lo,hi);return ix[e[ix]>start]
def describe(index,p,has_signal):
    s,e,prefix,v=index;point=candidates(index,p,p+1);left=max(0,p-100);right=p+101;around=candidates(index,left,right)
    if len(point):distance=0
    else:
        k=int(np.searchsorted(s,p,side='right'));poss=[]
        if k>0:poss.append(p-int(prefix[k-1])+1)
        if k<len(s):poss.append(int(s[k])-p)
        distance=min(poss) if poss else 1000000
    covered=0;end=-1
    for j in around:
        a,b=max(left,int(s[j])),min(right,int(e[j]));covered+=max(0,b-max(a,end));end=max(end,b)
    signal=lambda ix:float(np.log1p(v[ix].max())) if len(ix) else 0.
    return [1.,float(len(point)>0),signal(point) if has_signal else np.nan,float(np.log1p(min(distance,1000000))),covered/(right-left),signal(around) if has_signal else np.nan]

# Independent boundary/nesting fixture catches half-open interval errors.
fixture=(np.array([10,15,40]),np.array([20,25,50]),np.array([20,25,50]),np.array([1.,2.,3.]))
assert describe(fixture,20,True)[1]==1 and describe(fixture,25,True)[1]==0
assert abs(describe(fixture,25,True)[3]-np.log(2))<1e-12
assert abs(describe(fixture,30,True)[4]-25/131)<1e-12
sources=json.loads((OUT/'reference_chromatin_sources.json').read_text());soft=gzip.decompress((WORK/'GSE113480_family.soft.gz').read_bytes()).decode();assert 'Genome_build: hg19' in soft
to19=LiftOver(str(ROOT/'work/hg38ToHg19.over.chain.gz'));to38=LiftOver(str(ROOT/'work/hg19ToHg38.over.chain.gz'))
manifest=[]
for stem in ['salomon','siraj']:manifest.extend(csv.DictReader((WORK/(stem+'_manifest_local.csv')).open()))
coordinates=[];unmapped=0
for r in manifest:
    chrom='chr'+r['chromosome'];p=int(r['position_hg38'])-1
    if r['study']=='Siraj_OL42_SKNSH':
        c,q,_,_=r['variant'].split(':');old=(c,int(q)-1)
    else:
        mapped=to19.convert_coordinate(chrom,p);old=None
        if len(mapped)==1 and mapped[0][0]==chrom:
            c,q,strand,_=mapped[0];back=to38.convert_coordinate(c,q)
            if len(back)==1 and back[0][:3]==(chrom,p,strand):old=(c,q)
        if old is None:unmapped+=1
    coordinates.append({'hg38':(chrom,p),'hg19':old})
columns=[];matrices=[];checks=[]
for source in sources['files']:
    path=Path(source['local_file']);assert sha(path)==source['sha256'];label=source.get('label',source['accession']);hg19=source['accession']=='GSE113480';assembly='hg19' if hg19 else 'hg38';signal=source.get('file_type')!='footprints';index=load_bed(path,signal);x=np.full((len(manifest),6),np.nan,dtype='float32')
    for i,point in enumerate(coordinates):
        point=point[assembly]
        if point is None:x[i,0]=0;continue
        c,p=point
        if c in index:x[i]=describe(index[c],p,signal)
        else:x[i]=[1,0,0 if signal else np.nan,np.log1p(1000000),0,0 if signal else np.nan]
    # Independently verify exact point membership on sampled real coordinates.
    rng=np.random.default_rng(1209);tested=0
    for i in rng.choice(len(manifest),min(20,len(manifest)),replace=False):
        point=coordinates[i][assembly]
        if point is None or point[0] not in index:continue
        s,e,prefix,v=index[point[0]];actual=bool(np.any((s<=point[1]) & (point[1]<e)));assert bool(x[i,1])==actual;tested+=1
    columns.extend(label+'__'+name for name in ['coordinate_available','point_called_region','point_log1p_peak_signal','log1p_nearest_region_distance','called_region_fraction_201bp','window_log1p_max_peak_signal']);matrices.append(x);checks.append({'profile':label,'assembly':assembly,'source_sha256':sha(path),'real_point_membership_checks':tested,'point_called_regions':int(np.nansum(x[:,1]))})
features=np.concatenate(matrices,axis=1);target=WORK/'reference_chromatin_features_local.npz';np.savez_compressed(target,features=features)
receipt={'rows':len(manifest),'feature_count':features.shape[1],'columns':columns,'profiles':checks,'hg38_to_hg19_unresolved_Salomon':unmapped,'feature_sha256':sha(target),'manifest_sha256':{stem:sha(WORK/(stem+'_manifest_local.csv')) for stem in ['salomon','siraj']},'coordinate_convention':'Sources 0-based half-open; SNP input positions 1-based converted to 0-based. GSE hg19 queries use source Siraj coordinates or unique reciprocal Salomon liftover.','missingness':'Unavailable mapping is NaN plus coordinate_available=0. No called peak is a separate indicator; peak signal encoding zero outside called peaks does not mean measured zero accessibility. Footprint signal columns are unavailable, not interpreted as a quantitative measurement.','limitations':['Public bulk reference profiles; not single-cell evidence or the original Salomon donors','SK-N-SH is a neuroblastoma cell line; GSE profiles are induced/derived neural cells and fetal astrocytes, not aging patient measurements','ENCODE peaks and footprints derive from the same biological replicate; they are different covariates, not independent observations','These public experimental profiles can overlap foundation-model training resources; independence of AlphaGenome is not established','No reporter outcome-derived quantities or clinical traits enter these features'],'test_labels_used':False}
(OUT/'reference_chromatin_feature_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2),flush=True)
