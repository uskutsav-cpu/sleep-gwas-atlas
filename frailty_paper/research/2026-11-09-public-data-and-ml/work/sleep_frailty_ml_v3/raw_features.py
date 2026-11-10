"""Prediction features for new paired DNA, without reading assay labels.

The only AlphaGenome information accepted here is the audited local permissive
AVI score table. No measured reporter activity/effect enters these features.
"""
import sys,gzip,collections
from common import *
from inference_records import StoredFeatureInference
from motif_model import load_motifs,scan,mx
sys.path.insert(0,str(OLD));import model,model_mlx
sys.path.insert(0,str(ROOT/'work/stage5_env/lib/python3.13/site-packages'));from pyliftover import LiftOver

def load_bed(path,signal):
    grouped=collections.defaultdict(list)
    with gzip.open(path,'rt') as f:
        for line in f:
            if line.startswith(('#','track','browser')):continue
            row=line.rstrip().split('\t');s,e=int(row[1]),int(row[2]);assert 0<=s<e
            v=float(row[6]) if signal else np.nan;assert not signal or np.isfinite(v) and v>=0;grouped[row[0]].append((s,e,v))
    result={}
    for c,rows in grouped.items():
        a=np.asarray(sorted(rows),dtype='float64');s=a[:,0].astype('int64');e=a[:,1].astype('int64');result[c]=(s,e,np.maximum.accumulate(e),a[:,2])
    return result

def describe(index,p,signal):
    s,e,prefix,v=index
    def candidates(left,right):
        lo=np.searchsorted(prefix,left,side='right');hi=np.searchsorted(s,right,side='left');ix=np.arange(lo,hi);return ix[e[ix]>left]
    point=candidates(p,p+1);left=max(0,p-100);right=p+101;around=candidates(left,right)
    if len(point):distance=0
    else:
        k=int(np.searchsorted(s,p,side='right'));poss=[]
        if k>0:poss.append(p-int(prefix[k-1])+1)
        if k<len(s):poss.append(int(s[k])-p)
        distance=min(poss) if poss else 1000000
    covered=0;end=-1
    for j in around:
        a,b=max(left,int(s[j])),min(right,int(e[j]));covered+=max(0,b-max(a,end));end=max(end,b)
    sig=lambda ix:float(np.log1p(v[ix].max())) if len(ix) else 0.
    return [1.,float(len(point)>0),sig(point) if signal else np.nan,float(np.log1p(min(distance,1000000))),covered/(right-left),sig(around) if signal else np.nan]

class PublicReferenceChromatin:
    def __init__(self):
        sources=json.loads((OUT/'reference_chromatin_sources.json').read_text());self.profiles=[]
        for r in sources['files']:
            p=WORK/'reference_chromatin'/Path(r['local_file']).name
            if not p.is_file():p=Path(r['local_file'])
            assert sha(p)==r['sha256'];signal=r.get('file_type')!='footprints';self.profiles.append(('hg19' if r['accession']=='GSE113480' else 'hg38',signal,load_bed(p,signal)))
        self.to19=LiftOver(str(ROOT/'work/hg38ToHg19.over.chain.gz'));self.to38=LiftOver(str(ROOT/'work/hg19ToHg38.over.chain.gz'))
    def annotate(self,records):
        output=[]
        for r in records:
            c='chr'+str(int(str(r['chromosome']).removeprefix('chr')));assert 1<=int(c[3:])<=22;p=int(r['position_hg38'])-1;assert p>=0
            old=None;m=self.to19.convert_coordinate(c,p)
            if len(m)==1 and m[0][0]==c:
                cc,q,strand,_=m[0];back=self.to38.convert_coordinate(cc,q)
                if len(back)==1 and back[0][:3]==(c,p,strand):old=(cc,q)
            coords={'hg38':(c,p),'hg19':old};features=[]
            for assembly,signal,index in self.profiles:
                point=coords[assembly]
                if point is None:values=[0.,np.nan,np.nan,np.nan,np.nan,np.nan]
                elif point[0] in index:values=describe(index[point[0]],point[1],signal)
                else:values=[1.,0.,0. if signal else np.nan,np.log1p(1000000),0.,0. if signal else np.nan]
                features.extend(values)
            output.append(features)
        return np.asarray(output,dtype='float32')

def dna_features(a,b):
    def counts(x):
        x=x.astype('int32')-7;out=[]
        for k in [1,2,3]:
            ids=sum(x[:,o:x.shape[1]-k+1+o]*(4**(k-o-1)) for o in range(k));z=np.zeros((len(x),4**k),dtype='float32');np.add.at(z,(np.arange(len(x))[:,None],ids),1);out.append(z)
        return np.concatenate(out,1)
    def one(x,y):
        p,q=counts(x),counts(y);den=np.concatenate([np.full(4**k,x.shape[1]-k+1) for k in [1,2,3]]);return np.c_[(p+q)/2/den,np.abs(p-q)]
    return (.5*(one(a,b)+one(17-a[:,::-1],17-b[:,::-1]))).astype('float32')

class RawContrastFeatures:
    def __init__(self,motif,hyena,chrom):self.M=motif.shape[1]//3;self.motif=motif;self.hyena=hyena;self.chromatin=chrom
    def load_transform(self,z):self.mean=z['context_mean'];self.sd=z['context_sd'];self.vscale=z['variant_scale']
    def transform(self,indices):
        m=self.motif[indices];h=self.hyena[indices];c=np.c_[m[:,:self.M],m[:,2*self.M:],h[:,:128],self.chromatin[indices]];v=np.c_[m[:,self.M:2*self.M],h[:,128:]]
        return np.clip((np.where(np.isfinite(c),c,self.mean)-self.mean)/self.sd,-20,20).astype('float32'),np.clip(v/self.vscale,-20,20).astype('float32')

class RawDNAInference(StoredFeatureInference):
    def __init__(self,records,reference=None):
        assert records;self.records=records;self.cache={};self.manifest=[];self.data={'study':[]};refs=[];alts=[];positions=[]
        for r in records:
            a=r['ref_sequence'].upper();b=r['alt_sequence'].upper();assert len(a)==len(b) and len(a) in [200,270],'Use a200- or270-base paired construct, matching development data';assert set(a+b)<=set('ACGT')
            p=[j for j in range(len(a)) if a[j]!=b[j]];assert len(p)==1,'Exactly one substituted DNA base is required';p=p[0]
            c=str(int(str(r['chromosome']).removeprefix('chr')));assert 1<=int(c)<=22 and int(r['position_hg38'])>=1
            strand=r.get('genomic_strand','+');assert strand in ['+','-'];complement=str.maketrans('ACGT','TGCA');ref=a[p] if strand=='+' else a[p].translate(complement);alt=b[p] if strand=='+' else b[p].translate(complement)
            # Coordinates describe this substitution; engineered flanks need not be native genome.
            key=f'chr{c}:{int(r["position_hg38"])}:{ref}:{alt}'
            if 'variant_GRCh38' in r:assert r['variant_GRCh38']==key
            assert r['assay_context'] in STUDIES;self.data['study'].append(STUDIES.index(r['assay_context']));self.manifest.append({'variant_GRCh38':key});refs.append(np.array(['ACGT'.index(x)+7 for x in a],dtype='int32'));alts.append(np.array(['ACGT'.index(x)+7 for x in b],dtype='int32'));positions.append(p)
        self.data['study']=np.asarray(self.data['study'],dtype='int32');self.avi=None
        receipt=json.loads((OUT/'feature_receipt.json').read_text());assert sha(OLD/'hyenadna/model.safetensors')==receipt['datasets'][0]['base_weights_sha256']
        for r in json.loads((OUT/'motif_source_receipt.json').read_text())['sources']:assert sha(WORK/'motifs'/r['file'])==r['sha256']
        mx.metal.set_memory_limit(2*1024**3);mx.metal.set_cache_limit(64*1024**2)
        weights={k:mx.array(v) for k,v in model.read_weights().items()};groups,labels=load_motifs();groups={k:mx.array(v) for k,v in groups.items()};hyena=np.empty((len(records),512),dtype='float32');motif=np.empty((len(records),3*len(labels)),dtype='float32');dna=np.empty((len(records),168),dtype='float32')
        for length in [200,270]:
            indices=np.array([i for i,a in enumerate(refs) if len(a)==length],dtype='int64')
            for start in range(0,len(indices),16):
                ix=indices[start:start+16];a=np.stack([refs[i] for i in ix]);b=np.stack([alts[i] for i in ix]);p=np.array([positions[i] for i in ix]);hyena[ix]=np.array(model_mlx.strand_features(mx.array(a),mx.array(b),mx.array(p),weights));motif[ix]=np.array(scan(mx.array(a),mx.array(b),mx.array(p),groups)).astype('float16').astype('float32');dna[ix]=dna_features(a,b)
        reference=reference or PublicReferenceChromatin();chrom=reference.annotate(records);self.features=RawContrastFeatures(motif,hyena,chrom);self.original={'hyena_original':hyena,'dna_pattern':dna}
    def motif(self,indices):return self.features.motif[indices]
    def avi_features(self):
        if self.avi is None:
            lookup={};self.avi_source_receipts=[]
            # The v2 permissive lookup also contains the33 project variants,
            # which are not necessarily present in the v3 training manifests.
            # Never modify the immutable v3 training-target snapshot to add them.
            for directory,stem in [(OUT,'avi'),(ROOT/'outputs/sleep_frailty_ai_v2','avi')]:
                path=directory/(stem+'_exact_allele_scores.csv');audit=json.loads((directory/(stem+'_retrieval_audit.json')).read_text());assert not audit['failed_groups'] and audit['output_sha256']==sha(path);self.avi_source_receipts.append({'file':str(path),'sha256':sha(path)})
                for row in csv.DictReader(path.open()):
                    if row['status']!='exact_allele_match':continue
                    key=row['variant_GRCh38']
                    if key in lookup:assert all(abs(float(lookup[key][k])-float(row[k]))<1e-12 for k in ['avi_raw','avi_phred']),(key,'Inconsistent AVI cache values')
                    lookup[key]=row
            values=[]
            for r in self.manifest:
                row=lookup.get(r['variant_GRCh38']);values.append([float(row['avi_raw']),float(row['avi_phred']),1.] if row is not None and row['status']=='exact_allele_match' else [np.nan,np.nan,0.])
            self.avi=np.asarray(values,dtype='float32')
        return self.avi

