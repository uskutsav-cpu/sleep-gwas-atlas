"""Only encode added actual DNA; reuse all preserved old features read-only."""
import sqlite3
from v9_common import *
from legnet_runtime import base_model, encode
from extract_variant_local_features import local_pool
import torch


def main():
    assert_previous_sources()
    audit = json.loads((OUT/'revealed_source_augmentation_audit.json').read_text())
    rows = list(csv_rows(DATA/'extra_sequence_bank.csv.gz')); n = len(rows); b0 = audit['old_feature_bank_offset']
    lut = np.full(256, 255, dtype='uint8'); lut[[ord(c) for c in 'AGCT']] = [0,1,2,3]
    ref = np.empty((n,150), dtype='uint8'); alt = ref.copy(); sequences = []
    for i,r in enumerate(rows):
        assert int(r['bank_id']) == b0+i
        a,b = r['reference'],r['alternate']; p = int(r['changed_position']); start = max(0,min(p-75,len(a)-150))
        a,b = a[start:start+150],b[start:start+150]; p -= start
        assert len(a)==len(b)==150 and [j for j,(x,y) in enumerate(zip(a,b)) if x!=y]==[p]
        ref[i]=lut[np.frombuffer(a.encode(),dtype='uint8')];alt[i]=lut[np.frombuffer(b.encode(),dtype='uint8')];sequences.append((a,b,p))
    np.save(DATA/'extra_tokens_reference150.npy',ref);np.save(DATA/'extra_tokens_alternate150.npy',alt)
    device='mps' if torch.backends.mps.is_available() else 'cpu'; model=base_model().to(device)
    g=np.empty((n,2,128),dtype='float32');l=g.copy()
    for start in range(0,n,256):
        chunk=sequences[start:start+256];a=[x[0] for x in chunk];b=[x[1] for x in chunk];pp=[x[2] for x in chunk];x=encode(a+b,device);p=torch.tensor(pp+pp,device=device)
        with torch.no_grad():
            z=model.main(model.stem(x));gg=z.mean(-1).cpu().numpy();ll=local_pool(z,p).cpu().numpy()
        k=len(chunk);g[start:start+k,0]=gg[:k];g[start:start+k,1]=gg[k:];l[start:start+k,0]=ll[:k];l[start:start+k,1]=ll[k:]
    assert np.isfinite(g).all() and np.isfinite(l).all()
    np.save(DATA/'extra_global128.npy',g);np.save(DATA/'extra_local128.npy',l)
    db=sqlite3.connect((V7/'avi_exact_scores.sqlite').as_uri()+'?mode=ro',uri=True);av=[]
    for r in rows:
        s=db.execute('SELECT raw,phred,available FROM scores WHERE variant=?',(r['variant_GRCh38'],)).fetchone();assert s is not None,'Earlier exact source audit cache is missing this allele'
        av.append([s[0] if s[2] else np.nan,s[1] if s[2] else np.nan,s[2]])
    db.close();np.save(DATA/'extra_AVI.npy',np.array(av,dtype='float32'))
    outputs=['extra_tokens_reference150.npy','extra_tokens_alternate150.npy','extra_global128.npy','extra_local128.npy','extra_AVI.npy']
    write_json(OUT/'extra_feature_receipt.json',{'extra_physical_pairs':n,'original_bank_offset':b0,'outputs':{n:sha(DATA/n) for n in outputs},'original_feature_banks_copied_or_modified':False,'original_pretrained_checkpoint_sha256':sha(V8/'models/author_WTC11_original.safetensors'),'fresh_test_sequences_or_targets_used_for_training':False,'scope':'Frozen numerical DNA features for additional revealed experimental development rows; no successful transfer or disease-risk validation is implied.'})
    assert_previous_sources();print('ADDED DNA FEATURES READY; ORIGINAL BANKS UNCHANGED',n,flush=True)


if __name__=='__main__':main()
