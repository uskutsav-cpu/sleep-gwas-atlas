"""Prepare all eligible existing AVI rows plus novel native-context pilot rows.

Only predicted AVI targets are included here. Laboratory outcomes and actual
activity flags are deliberately absent. Novel teacher validation is separate
from legacy assay-development rows.
"""
import sys,time
from common import *
def main():
    retrieval=json.loads((OUT/'teacher_pilot_retrieval.json').read_text());assert not retrieval['failed'];regions=json.loads((DATA/'teacher_regions.json').read_text());selected=[r for r in regions if r['pilot']];old_audit=json.loads((V3/'avi_retrieval_audit.json').read_text());score_path=V3/'avi_exact_allele_scores.csv';assert old_audit['output_sha256']==sha(score_path);scores={r['variant_GRCh38']:r for r in csv.DictReader(score_path.open())};buckets={200:[],270:[]};legacy_receipts=[]
    for stem in ['salomon','siraj']:
        source=ROOT/'work/sleep_frailty_ml_v3'/(stem+'_sequences_local.npz');feat=source.with_name(stem+'_features_local.npz');receipt=next(r for r in json.loads((V3/'feature_receipt.json').read_text())['datasets'] if r['study']==stem);assert sha(source)==receipt['input_sha256'] and sha(feat)==receipt['features_sha256']
        with np.load(source) as z:ref=z['ref'];alt=z['alt'];pos=z['position'];chrom=z['chromosome']
        with np.load(feat) as z:features=z['hyena_original']
        manifest=list(csv.DictReader(source.with_name(stem+'_manifest_local.csv').open()));assert len(ref)==len(manifest)
        for i in np.flatnonzero(chrom>=4):
            s=scores[manifest[i]['variant_GRCh38']]
            if s['status']!='exact_allele_match':continue
            split=0 if chrom[i]>=6 else 2
            buckets[ref.shape[1]].append({'ref':ref[i],'alt':alt[i],'position':int(pos[i]),'target':[float(s['avi_raw']),float(s['avi_phred'])],'features':features[i],'split':split,'anchor_id':stem+'_'+manifest[i]['variant_GRCh38'],'source':stem,'chromosome':int(chrom[i])})
        legacy_receipts.append({'study':stem,'sequence_sha256':sha(source),'feature_sha256':sha(feat),'base_weights_sha256':receipt['base_weights_sha256']})
    new_rows=[];cached={r['anchor_id']:r for r in retrieval['chunks']}
    for r in selected:
        path=DATA/'teacher_cache'/(r['anchor_id']+'.json');assert sha(path)==cached[r['anchor_id']]['cache_sha256'];d=json.loads(path.read_text());assert d['complete'] and d['anchor']==r
        for p in d['pairs']:
            j=p['offset']+15;ref=d['native_sequence_300bp'][j:j+270];alt=ref[:134]+p['alternate_base']+ref[135:];assert len(ref)==270 and ref[134]==p['reference_base'] and sum(a!=b for a,b in zip(ref,alt))==1
            new_rows.append({'ref':np.array(['ACGT'.index(c)+7 for c in ref],dtype='uint8'),'alt':np.array(['ACGT'.index(c)+7 for c in alt],dtype='uint8'),'position':134,'target':[p['avi_raw'],p['avi_phred']],'split':0 if r['split']=='train' else 1,'anchor_id':r['anchor_id'],'source':'novel_native','chromosome':int(r['chromosome'][3:])})
    sys.path.insert(0,str(OLD/'runtime'));sys.path.insert(0,str(OLD));import mlx.core as mx;import model,model_mlx
    mx.metal.set_memory_limit(1536*1024**2);mx.metal.set_cache_limit(64*1024**2);weights={k:mx.array(v) for k,v in model.read_weights().items()};assert all(r['base_weights_sha256']==sha(OLD/'hyenadna/model.safetensors') for r in legacy_receipts);compiled=mx.compile(lambda a,b,p:model_mlx.strand_features(a,b,p,weights));feature_dir=DATA/'teacher_feature_shards';feature_dir.mkdir(exist_ok=True);signature=sha(DATA/'teacher_regions.json')+'_'+sha(OUT/'teacher_pilot_retrieval.json')+'_'+sha(OLD/'model_mlx.py')+'_'+sha(OLD/'hyenadna/model.safetensors');started=time.time()
    for start in range(0,len(new_rows),512):
        end=min(len(new_rows),start+512);path=feature_dir/f'{start}_{end}.npz'
        if path.exists():
            with np.load(path) as z:assert str(z['signature'])==signature;f=z['features']
        else:
            values=[]
            for j in range(start,end,32):
                rr=new_rows[j:min(end,j+32)];a=np.stack([r['ref'] for r in rr]).astype('int32');b=np.stack([r['alt'] for r in rr]).astype('int32');p=np.array([134]*len(rr),dtype='int32');values.append(np.array(compiled(mx.array(a),mx.array(b),mx.array(p))))
            f=np.concatenate(values);np.savez_compressed(path,signature=np.array(signature),features=f)
        assert f.shape==(end-start,512) and np.isfinite(f).all()
        for r,x in zip(new_rows[start:end],f):r['features']=x
        if start%4096==0:print('Novel teacher embeddings',end,'/',len(new_rows),'seconds',round(time.time()-started,1),flush=True)
    buckets[270].extend(new_rows);receipt=[]
    for length,rr in buckets.items():
        path=DATA/f'teacher_data_{length}.npz';arr={k:np.stack([r[k] for r in rr]) for k in ['ref','alt','position','target','features','split','anchor_id','source','chromosome']};arr['target']=arr['target'].astype('float32');arr['position']=arr['position'].astype('int32');arr['split']=arr['split'].astype('int32');arr['chromosome']=arr['chromosome'].astype('int32');assert np.all(np.isfinite(arr['target'])) and np.all(arr['target'][:,1]>=0);np.savez_compressed(path,**arr);receipt.append({'length':length,'file':str(path),'sha256':sha(path),'rows':len(rr),'training_rows':int(np.sum(arr['split']==0)),'novel_validation_rows':int(np.sum(arr['split']==1)),'legacy_development_probe_rows':int(np.sum(arr['split']==2))})
    write_json(OUT/'teacher_data_receipt.json',{'datasets':receipt,'legacy_sources':legacy_receipts,'new_exact_teacher_rows':len(new_rows),'legacy_AVI_sha256':sha(score_path),'base_weights_sha256':sha(OLD/'hyenadna/model.safetensors'),'feature_code_sha256':sha(OLD/'model_mlx.py'),'teacher_targets':'AVIraw and PHRED predictions only; no lab call/activity/frailty/insomnia labels in these files','split_codes':{'0':'Train chromosomes6–22','1':'Novel validation anchors on4–5; used for teacher checkpoint selection','2':'Legacy4–5assay-region AVI probe; previously used development data, not primary selection'},'reserved_teacher_test_loaded':False,'V3_model_unchanged':sha(V3/'final_model_registry.json')==json.loads((OUT/'teacher_region_plan.json').read_text())['old_V3_registry_sha256']});print('TEACHER DATA READY',json.dumps(receipt),flush=True)
if __name__=='__main__':main()
