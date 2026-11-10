"""Resumable new local features; verified cached context cannot encode labels."""
import gzip
from common import *
from sequence_features import Encoder,SLICES
def main():
    assert not (OUT/'model_search_started.json').exists()
    audit=json.loads((OUT/'development_data_audit.json').read_text())
    for name,h in audit['outputs'].items():assert sha(DATA/name)==h,name
    with gzip.open(DATA/'sequence_bank.csv.gz','rt') as f:bank=list(csv.DictReader(f))
    previous={'v6':json.loads((V6OUT/'fixed_features_receipt.json').read_text()),'kosicki':json.loads((V6OUT/'kosicki_fixed_features_receipt.json').read_text()),'hepg2':json.loads((V6OUT/'additional_study_fixed_features_receipt.json').read_text())}
    paths={'v6':Path(previous['v6']['feature_path']),'kosicki':V6/'kosicki_fixed_features.npy','hepg2':V6/'additional_study_fixed_features.npy'}
    receipts={};cached={}
    for kind,p in paths.items():
        h=previous[kind].get('feature_sha256',previous[kind].get('sha256'));assert sha(p)==h,p;receipts[kind]={'path':str(p),'sha256':h};cached[kind]=np.load(p,mmap_mode='r')
    signature=hashlib.sha256((sha(DATA/'sequence_bank.csv.gz')+sha(WORK/'sequence_features.py')+sha(OLD/'model_mlx.py')+sha(OLD/'hyenadna/model.safetensors')+sha(V6OUT/'selected_motifs.json')).encode()).hexdigest();path=DATA/'feature_bank.npy';progress=OUT/'feature_progress.json'
    if path.exists() and progress.exists():
        saved=json.loads(progress.read_text());assert saved['signature']==signature;start=saved['completed_rows'];x=np.load(path,mmap_mode='r+')
    else:
        start=0;x=np.lib.format.open_memmap(path,mode='w+',dtype='float32',shape=(len(bank),1145));write_json(progress,{'signature':signature,'completed_rows':0,'total':len(bank)})
    encoder=Encoder();started=time.time()
    for begin in range(start,len(bank),1024):
        end=min(begin+1024,len(bank));rows=bank[begin:end];contexts=np.zeros((len(rows),128),dtype='float32');unknown=[]
        for j,r in enumerate(rows):
            kind=r['context_feature_source']
            if kind in cached:contexts[j]=cached[kind][int(r['context_feature_index']),:128]
            else:unknown.append(j)
        value=encoder.features([r['reference'] for r in rows],[r['alternate'] for r in rows],[int(r['changed_position']) for r in rows],context=contexts)
        if unknown:
            rr=[rows[j] for j in unknown];full=encoder.features([r['reference'] for r in rr],[r['alternate'] for r in rr],[int(r['changed_position']) for r in rr]);value[unknown,:128]=full[:,:128]
        x[begin:end]=value;x.flush();write_json(progress,{'signature':signature,'completed_rows':end,'total':len(bank),'seconds_this_run':time.time()-started})
        if begin%8192==0 or end==len(bank):print('V7 measured-effect features',end,'/',len(bank),'seconds',round(time.time()-started,1),flush=True)
    write_json(OUT/'feature_receipt.json',{'signature':signature,'path':str(path),'sha256':sha(path),'rows':len(bank),'columns':1145,'slices':SLICES,'source_context_caches':receipts,'backbone_sha256':sha(OLD/'hyenadna/model.safetensors'),'motif_selection_sha256':sha(V6OUT/'selected_motifs.json'),'feature_code_sha256':sha(WORK/'sequence_features.py'),'targets_or_AVI_used_in_sequence_features':False,'scope':'Original fixed DNA backbone. New variant51bp evidence, literal51bpsequence counts and51bp motifs. Context is actual native sequence up to151bp; old byte-verified128dimensional context means reused. No measured expression/activity values or taught V4 weights are inputs.'});print('V7 SEQUENCE FEATURES COMPLETE',len(bank),flush=True)
if __name__=='__main__':main()
