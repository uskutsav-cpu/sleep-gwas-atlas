"""Reconstruct a separate measured neuronal assay for local research evaluation.

Source article CC BY-NC-ND 4.0: retain inputs and extracts locally; this script
does not redistribute sequence-level derivatives or model weights learned from
them. No external labels enter v2 training or selection.
"""
from pathlib import Path
import sys,re,json,csv,hashlib,collections
import numpy as np
from read_xlsx import dictionaries,rows
ROOT=Path(__file__).resolve().parents[2];WORK=Path(__file__).resolve().parent;OUT=ROOT/'outputs/sleep_frailty_ai_v2';SOURCE=ROOT/'work/retallick_neural/published_supp2.xlsx'
seq={};conflict=set()
for r in dictionaries(SOURCE,'SI Data 2.9'):
    name,s=r['SNP'],r['Sequence'].upper()
    if name in seq and seq[name]!=s:conflict.add(name)
    else:seq[name]=s
for k in conflict:seq.pop(k,None)
groups=collections.defaultdict(list)
for r in dictionaries(SOURCE,'SI_Data 2.12'):
    if r['type'] in ['ref','alt']:groups[r['ID']].append(r)
excluded=collections.Counter();out=[];ref=[];alt=[];ys=[];positions=[]
lut=np.zeros(256,dtype='uint8');lut[[ord(x) for x in 'ACGT']]=[7,8,9,10]
for ident,rs in groups.items():
    if len(rs)!=2 or {r['type'] for r in rs}!={'ref','alt'}:excluded['not_unique_pair']+=1;continue
    labels={r['emVar.group'] for r in rs}
    if len(labels)!=1:excluded['label_disagreement']+=1;continue
    label=next(iter(labels))
    if label not in ['fdr <= 0.05','n.s']:excluded['excluded_label']+=1;continue
    a,b=(next(r for r in rs if r['type']==t) for t in ['ref','alt']);m=re.fullmatch(r'chr(\d+):rs\d+:(\d+)',ident)
    if not m or not 1<=int(m[1])<=22:excluded['non_autosomal_id']+=1;continue
    sa,sb=seq.get(a['CRS']),seq.get(b['CRS'])
    if not sa or not sb or len(sa)!=199 or len(sb)!=199 or set(sa+sb)-set('ACGT'):excluded['sequence_unusable']+=1;continue
    change=[i for i,(x,y) in enumerate(zip(sa,sb)) if x!=y]
    if len(change)!=1:excluded['not_single_base_pair']+=1;continue
    p=change[0]
    if sa[p]!=a['allele'] or sb[p]!=b['allele']:excluded['source_allele_conflict']+=1;continue
    ref.append(lut[np.frombuffer(sa.encode(),dtype='uint8')]);alt.append(lut[np.frombuffer(sb.encode(),dtype='uint8')]);ys.append(int(label=='fdr <= 0.05'));positions.append(p)
    out.append({'variant':ident,'rsid':a['rsID'],'chromosome':int(m[1]),'source_position_build_unresolved':int(m[2]),'source_ref':sa[p],'source_alt':sb[p],'changed_position':p,'label':ys[-1]})
prior=json.loads((ROOT/'outputs/stage6_neuronal_mpra_kmer_baseline.json').read_text())
assert len(out)==prior['train_n']+prior['test_n'] and sum(ys)==prior['train_detected']+prior['test_detected'],(len(out),sum(ys))
np.savez_compressed(WORK/'external_retallick_local.npz',ref=np.stack(ref),alt=np.stack(alt),label=np.array(ys,dtype='float32'),position=np.array(positions,dtype='int32'))
with (WORK/'external_retallick_manifest_local.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(out[0]));w.writeheader();w.writerows(out)
receipt={'source':'Retallick-Townsley et al. 2026 Supplementary Data 2 baseline','doi':'10.1038/s41467-026-77114-x','source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),'eligible_pairs':len(out),'detected_calls':sum(ys),'source_groups':len(groups),'excluded':dict(excluded),'sequence_length':199,'scope':'External measured assay relative to Salomon fitting; no external labels enter model training or selection','current_limitations':['Genome build and exact-allele AVI mapping unresolved; must verify before annotation','199-base constructs differ from 270-base training contexts','Overlap with training regions must be checked before performance reporting','The broader project previously inspected aggregate results and a separate baseline on this dataset','AlphaGenome/base-model exposure to this assay is unresolved; not proven model-independent','Source article CC BY-NC-ND 4.0; do not redistribute sequence-level extracts'],'status':'reconstructed locally; mapping and leakage audit pending'}
(OUT/'external_retallick_preparation_audit.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2),flush=True)
