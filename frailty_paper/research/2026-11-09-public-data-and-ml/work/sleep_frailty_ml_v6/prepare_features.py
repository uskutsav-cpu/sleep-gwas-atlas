"""Fixed DNA encodings and declared motif features, with replayable shards."""
import sys,gzip,collections
from common import *
sys.path.insert(0,str(OLD/'runtime'));sys.path.insert(0,str(OLD));import model,model_mlx,mlx.core as mx
sys.path.insert(0,str(ROOT/'work/sleep_frailty_ml_v3'));import motif_model
TF_LIST='ARNTL CLOCK NPAS2 RORA RORB RORC NR1D1 NR1D2 BHLHE40 BHLHE41 PER1 CREB1 CREM ATF1 ATF2 ATF3 ATF4 ATF6 JUN JUNB JUND FOS FOSB FOSL1 FOSL2 MEF2A MEF2B MEF2C MEF2D MYOD1 MYOG MYF5 MYF6 PAX3 PAX7 SPI1 CEBPA CEBPB CEBPD CEBPE NFKB1 NFKB2 REL RELA RELB STAT1 STAT2 STAT3 STAT4 STAT5A STAT5B STAT6 NR3C1 ESR1 ESR2 AR CTCF FOXO1 FOXO3 FOXO4 HIF1A HSF1 REST NEUROD1 NEUROD2 NEUROG1 NEUROG2 ASCL1 ASCL2 TCF4 SOX2 SOX3 SOX6 SOX10 SOX11 POU3F2 POU3F3 POU3F4 RFX1 RFX2 RFX3 RFX4 RFX5 RFX7 GATA1 GATA2 GATA3 GATA4 GATA6 TAL1 RUNX1 RUNX2 RUNX3 IKZF1 IKZF2 IRF1 IRF2 IRF3 IRF4 IRF5 IRF7 IRF8 IRF9 HNF1A HNF4A HNF4G FOXA1 FOXA2 FOXA3 SREBF1 SREBF2 PPARG PPARA PPARD RXRA RXRB RXRG NFIA NFIB NFIC NFIX SP1 SP2 SP3 SP4 KLF1 KLF2 KLF3 KLF4 KLF5 KLF6 KLF9 KLF10 KLF13 KLF15 KLF16 EGR1 EGR2 EGR3 EGR4 ELK1 ELK3 ELK4 ETS1 ETS2 ETV1 ETV4 ETV5 ETV6 YY1 NRF1 NFE2L1 NFE2L2 NFE2 TP53 TP63 TP73 TBP'.split()
def selected_motifs():
    groups,labels=motif_model.load_motifs();pwm_by_name={};j=0
    for width,items in groups.items():
        for pwm in items:pwm_by_name[labels[j]['motif']]=pwm;j+=1
    chosen=[]
    for tf in TF_LIST:
        opts=[r for r in labels if r['human_gene']==tf]
        if opts:chosen.append(min(opts,key=lambda r:(r['quality'],r['length'],r['motif'])))
    chosen.sort(key=lambda r:(r['length'],r['motif']));out=collections.defaultdict(list)
    for r in chosen:out[r['length']].append(pwm_by_name[r['motif']])
    return {w:np.stack(v) for w,v in out.items()},chosen
def pattern(ref,alt):
    def counts(a):
        a=a.astype('int32')-7;parts=[]
        for k in [1,2,3]:
            ids=sum(a[:,o:a.shape[1]-k+1+o]*(4**(k-o-1)) for o in range(k));out=np.zeros((len(a),4**k),dtype='float32');np.add.at(out,(np.arange(len(a))[:,None],ids),1);parts.append(out)
        return np.concatenate(parts,1)
    def one(a,b):
        p,q=counts(a),counts(b);den=np.concatenate([np.full(4**k,a.shape[1]-k+1) for k in [1,2,3]])
        return np.c_[(p+q)/2/den,np.abs(p-q)]
    return (.5*(one(ref,alt)+one(17-ref[:,::-1],17-alt[:,::-1]))).astype('float32')
def main():
    if (DATA/'fixed_feature_bank.npy').is_symlink():raise RuntimeError('V6 shares immutable verified V5 encodings; rebuild into a new independent path rather than overwriting the shared cache.')
    assert not (OUT/'model_search_started.json').exists()
    audit=json.loads((OUT/'measured_data_audit.json').read_text())
    for name,h in audit['outputs'].items():assert sha(DATA/name)==h
    with gzip.open(DATA/'sequence_bank.csv.gz','rt') as f:bank=list(csv.DictReader(f))
    with np.load(DATA/'local151_sequences.npz') as z:ref=z['ref'];alt=z['alt'];pos=z['position']
    groups,labels=selected_motifs();write_json(OUT/'selected_motifs.json',{'requested_TFs':TF_LIST,'selection':'One best available quality then shortest human HOCOMOCO motif per named TF, chosen before fitting without outcomes','motifs':labels,'missing_TFs':sorted(set(TF_LIST)-{r['human_gene'] for r in labels}),'source_sha256':sha(ROOT/'work/sleep_frailty_ml_v3/motifs/H13CORE_pwms.txt')})
    dims={'local_hyena':[0,512],'full_hyena':[512,1024],'local_pattern':[1024,1192],'local_motif':[1192,1192+3*len(labels)]};D=dims['local_motif'][1]
    groups={k:mx.array(v) for k,v in groups.items()}
    weights={k:mx.array(v) for k,v in model.read_weights().items()};mx.metal.set_memory_limit(1536*1024**2);mx.metal.set_cache_limit(64*1024**2)
    encode=mx.compile(lambda a,b,p:model_mlx.strand_features(a,b,p,weights));scan=mx.compile(lambda a,b,p:motif_model.scan(a,b,p,groups))
    signature=hashlib.sha256((sha(DATA/'local151_sequences.npz')+sha(DATA/'sequence_bank.csv.gz')+sha(OLD/'hyenadna/model.safetensors')+sha(OUT/'selected_motifs.json')+sha(Path(__file__))).encode()).hexdigest()
    matrix=DATA/'fixed_feature_bank.npy';x=np.lib.format.open_memmap(matrix,mode='w+',dtype='float32',shape=(len(bank),D));shards=DATA/'feature_shards';shards.mkdir(exist_ok=True)
    started=time.time()
    for start in range(0,len(bank),1024):
        end=min(len(bank),start+1024);path=shards/f'local_{start}_{end}.npz'
        if path.exists():
            with np.load(path) as z:assert str(z['signature'])==signature;value=z['features']
        else:
            pieces=[]
            for j in range(start,end,32):
                k=min(j+32,end);a=mx.array(ref[j:k].astype('int32'));b=mx.array(alt[j:k].astype('int32'));p=mx.array(pos[j:k]);h=np.array(encode(a,b,p));m=np.array(scan(a,b,p));pieces.append(np.c_[h,pattern(ref[j:k],alt[j:k]),m])
            value=np.concatenate(pieces);assert np.isfinite(value).all();np.savez_compressed(path,signature=np.array(signature),features=value)
        x[start:end,:512]=value[:,:512];x[start:end,1024:]=value[:,512:]
        if start%8192==0:print('V5 fixed local encodings',end,'/',len(bank),'seconds',round(time.time()-started,1),flush=True)
    receipts=json.loads((ROOT/'outputs/sleep_frailty_ai_v3/feature_receipt.json').read_text())['datasets']
    for stem in ['salomon','siraj']:
        path=ROOT/'work/sleep_frailty_ml_v3'/(stem+'_features_local.npz');r=next(r for r in receipts if r['study']==stem);assert sha(path)==r['features_sha256']
        with np.load(path) as z:full=z['hyena_original']
        ix=np.array([int(r['bank_id']) for r in bank if r['sequence_source']==stem]);old=np.array([int(bank[i]['original_row_index']) for i in ix]);x[ix,512:1024]=full[old];del full
    by_length=collections.defaultdict(list)
    for r in bank:
        if r['sequence_source']=='kircher':by_length[len(r['reference'])].append(int(r['bank_id']))
    lut=np.zeros(256,dtype='uint8');lut[[ord(c) for c in 'ACGT']]=[7,8,9,10]
    for length,ids in by_length.items():
        for start in range(0,len(ids),32):
            ix=ids[start:start+32];a=np.stack([lut[np.frombuffer(bank[i]['reference'].encode(),dtype='uint8')] for i in ix]);b=np.stack([lut[np.frombuffer(bank[i]['alternate'].encode(),dtype='uint8')] for i in ix]);p=np.array([int(bank[i]['changed_position']) for i in ix],dtype='int32');x[ix,512:1024]=np.array(encode(mx.array(a.astype('int32')),mx.array(b.astype('int32')),mx.array(p)))
        print('Kircher full-context encoding length',length,'pairs',len(ids),flush=True)
    x.flush()
    for start in range(0,len(bank),4096):assert np.isfinite(x[start:start+4096]).all()
    write_json(OUT/'fixed_features_receipt.json',{'bank_rows':len(bank),'feature_count':D,'slices':dims,'feature_path':str(matrix),'feature_sha256':sha(matrix),'signature':signature,'source_base_sha256':sha(OLD/'hyenadna/model.safetensors'),'motif_count':len(labels),'source_code_sha256':sha(Path(__file__)),'targets_or_AlphaGenome_outputs_used':False,'scope':'Fixed original reference-pretrained DNA encoder, literal local sequence counts and declared DNA-binding motifs. No V3/V4 task-trained weights.','seconds':time.time()-started});print('V5 FIXED FEATURES COMPLETE',len(bank),D,flush=True)
if __name__=='__main__':main()
