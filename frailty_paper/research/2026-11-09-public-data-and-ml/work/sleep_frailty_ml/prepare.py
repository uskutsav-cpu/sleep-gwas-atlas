from pathlib import Path
import csv,gzip,json,hashlib,collections
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'outputs/sleep_frailty_ai_v1'
WORK=Path(__file__).resolve().parent
OUT.mkdir(exist_ok=True)
SOURCE=ROOT/'work/salomon_neuronal_variant_readout.tsv.gz'
comp=str.maketrans('ACGT','TGCA')
def rc(s):return s.translate(comp)[::-1]
def canonical(s):return min(s,rc(s))

def prepare():
    with gzip.open(SOURCE,'rt') as f:rows=list(csv.DictReader(f,delimiter='\t'))
    counts=collections.Counter(r['chr_pos_ref_alt'] for r in rows)
    kept=[];excluded=collections.Counter()
    for r in rows:
        key=r['chr_pos_ref_alt']; parts=key.split('-')
        if counts[key]!=1:excluded['duplicated_genomic_contrast']+=1;continue
        if len(parts)!=4 or not parts[0].isdigit() or not 1<=int(parts[0])<=22:excluded['nonautosomal_or_malformed']+=1;continue
        if r['NGN2_variant_has_readout']!='True' or r['NGN2_variant_is_significant'] not in ['True','False']:excluded['readout_or_label_unavailable']+=1;continue
        ref,alt=r['ref_seq'].upper(),r['alt_seq'].upper()
        if len(ref)!=270 or len(alt)!=270 or set(ref+alt)-set('ACGT'):excluded['unusable_sequence']+=1;continue
        changes=[i for i,(a,b) in enumerate(zip(ref,alt)) if a!=b]
        if len(changes)!=1:excluded['not_single_substitution']+=1;continue
        p=changes[0]; a,b=parts[2:]
        if len(a)!=1 or len(b)!=1 or r['strand'] not in ['+','-']:excluded['allele_or_strand_unusable']+=1;continue
        if r['strand']=='-':a,b=rc(a),rc(b)
        # Input variant_pos is a zero-based genomic offset, not always sequence offset.
        expected_p=int(parts[1])-1-int(r['start'])
        if r['strand']=='-':expected_p=269-expected_p
        if p!=expected_p or ref[p]!=a or alt[p]!=b:excluded['coordinate_or_allele_conflict']+=1;continue
        chrom=int(parts[0]);split='train' if chrom<=16 else ('validation' if chrom<=18 else 'test')
        kept.append({'variant':key,'chromosome':chrom,'split':split,'reference':ref,'alternate':alt,'changed_position':p,'label':int(r['NGN2_variant_is_significant']=='True'),'gene':r['associated_gene_name']})
    # Reference-sequence matching also catches same enhancer with different variants.
    seen={s:set() for s in ['train','validation','test']}
    for r in kept:seen[r['split']].update([canonical(r['reference']),canonical(r['alternate'])])
    filtered=[]
    for r in kept:
        later=seen['test'] | (seen['validation'] if r['split']=='train' else set())
        if r['split']!='test' and any(canonical(r[k]) in later for k in ['reference','alternate']):excluded['exact_sequence_overlap_with_later_split']+=1;continue
        filtered.append(r)
    kept=filtered
    lut=np.zeros(256,dtype=np.uint8);lut[[ord(x) for x in 'ACGT']]=[7,8,9,10]
    def tok(s):return lut[np.frombuffer(s.encode(),dtype=np.uint8)]
    ref=np.stack([tok(r['reference']) for r in kept]);alt=np.stack([tok(r['alternate']) for r in kept])
    data={'ref':ref,'alt':alt,'position':np.array([r['changed_position'] for r in kept],dtype=np.int32),'label':np.array([r['label'] for r in kept],dtype=np.float32),'chromosome':np.array([r['chromosome'] for r in kept],dtype=np.int32)}
    np.savez_compressed(WORK/'prepared.npz',**data)
    # Sequence-level data are retained locally only. Publish source identifiers and splits.
    with (OUT/'split_manifest.csv').open('w',newline='') as f:
        fields=['variant','chromosome','split','changed_position','label','gene'];w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows([{k:r[k] for k in fields} for r in kept])
    summary={'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),'source_rows':len(rows),'excluded':dict(excluded),'eligible':len(kept),'splits':{s:{'n':sum(r['split']==s for r in kept),'detected_calls':sum(r['label'] for r in kept if r['split']==s),'chromosomes':sorted({r['chromosome'] for r in kept if r['split']==s})} for s in seen},'label':'Source detected allelic reporter effect versus nonsignificant; not biological inactive/active truth or disease status','input_features':'Only paired DNA; no AlphaGenome, Enformer, outcome-derived annotations or GWAS p values'}
    (OUT/'data_audit.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2),flush=True)
    return data

if __name__=='__main__':prepare()
