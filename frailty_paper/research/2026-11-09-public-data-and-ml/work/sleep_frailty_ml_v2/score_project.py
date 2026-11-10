"""Score the 33 public project contrasts with the validation-selected model."""
from pathlib import Path
import sys,csv,json
import numpy as np
from common import *
from inference import infer_record
from prepare_features import kmer_features
rows=list(csv.DictReader((OLD/'project_sequence_pairs.csv').open()));assert len(rows)==33
avi={r['variant_GRCh38']:r for r in csv.DictReader((OUT/'avi_exact_allele_scores.csv').open())};aa=[]
for r in rows:
    chrom=r['chromosome'];chrom=chrom if chrom.startswith('chr') else 'chr'+chrom
    key=f"{chrom}:{r['position_hg38']}:{r['genome_reference_allele']}:{r['tested_alternate_allele']}";a=avi[key]
    aa.append([float(a['avi_raw']),float(a['avi_phred']),1] if a['status']=='exact_allele_match' else [np.nan,np.nan,0])
aa=np.array(aa,dtype='float32');lut=np.zeros(256,dtype='int32');lut[[ord(c) for c in 'ACGT']]=[7,8,9,10]
ref=np.stack([lut[np.frombuffer(r['reference_sequence'].encode(),dtype='uint8')] for r in rows]);alt=np.stack([lut[np.frombuffer(r['alternate_sequence'].encode(),dtype='uint8')] for r in rows]);pos=np.array([np.flatnonzero(a!=b)[0] for a,b in zip(ref,alt)],dtype='int32')
sys.path.insert(0,str(OLD/'runtime'));sys.path.insert(0,str(OLD));import mlx.core as mx;import model,model_mlx
w={k:mx.array(v) for k,v in model.read_weights().items()};emb=np.array(model_mlx.strand_features(mx.array(ref),mx.array(alt),mx.array(pos),w));km=kmer_features(ref,alt)
features={'avi':aa,'kmer':km,'kmer_avi':np.c_[km,aa],'hyena':emb,'hyena_avi':np.c_[emb,aa]};records={r['name']:r for r in json.loads((OUT/'search_results.json').read_text())};selection=json.loads((OUT/'model_selection.json').read_text());selected=selection['selected_model'];pred=infer_record(records[selected],features,records,rows)
target=OUT/'project_selected_model_predictions.csv'
with target.open('w',newline='') as f:
    fields=['variant','selected_model','uncalibrated_assay_call_score','avi_raw','avi_phred','interpretation'];writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader()
    for r,p,a in zip(rows,pred,aa):writer.writerow({'variant':r['variant'],'selected_model':selected,'uncalibrated_assay_call_score':float(p),'avi_raw':float(a[0]) if np.isfinite(a[0]) else '', 'avi_phred':float(a[1]) if np.isfinite(a[1]) else '','interpretation':'Exploratory neuronal reporter ranking; no measured candidate label, clinical risk, disease association or causal-gene proof'})
(OUT/'project_prediction_receipt.json').write_text(json.dumps({'rows':len(rows),'selected_model':selected,'selection_sha256':sha(OUT/'model_selection.json'),'input_sha256':sha(OLD/'project_sequence_pairs.csv'),'output_sha256':sha(target),'candidate_accuracy_known':False},indent=2)+'\n');print('Scored project contrasts',len(rows),flush=True)
