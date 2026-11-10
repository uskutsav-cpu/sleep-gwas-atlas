"""Prepare additional measured Siraj constructs and strict cross-study splits.

Does not predict on the held-out Siraj test set or fit a model.
"""
from pathlib import Path
import sys,csv,gzip,zipfile,json,hashlib,collections,math,bisect
import numpy as np
ROOT=Path(__file__).resolve().parents[2];WORK=Path(__file__).resolve().parent;OUT=ROOT/'outputs/sleep_frailty_ai_v3';OLD=ROOT/'work/sleep_frailty_ml'
sys.path.insert(0,str(ROOT/'work/stage5_env/lib/python3.13/site-packages'))
from pyliftover import LiftOver
SOURCE=Path('/Volumes/Extreme SSD/codex-sleep-frailty-work/siraj2026_neural/OL42_SKNSH_emVAR_glm_paired_20221220.out.gz');TABLE=Path('/Volumes/Extreme SSD/codex-sleep-frailty-work/siraj2026_tables.zip')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert sha(SOURCE)=='cdb42e4b7169dea14ca3db543a83764a0b849ad2acdd8aee8734af63caf5d282'
out37=LiftOver(str(ROOT/'work/hg19ToHg38.over.chain.gz'));back=LiftOver(str(ROOT/'work/hg38ToHg19.over.chain.gz'));comp=str.maketrans('ACGT','TGCA')
def rc(s):return s.translate(comp)[::-1]
def canonical(s):return min(s,rc(s))
def split(chrom):return 'test' if chrom<=3 else 'validation' if chrom<=5 else 'train'
with gzip.open(SOURCE,'rt') as f:raw=[r for r in csv.DictReader(f,delimiter='\t') if r['comb'].endswith('_center_fwd_ref')]
counts=collections.Counter(r['SNP'] for r in raw);wanted={'chr'+r['SNP'] for r in raw};sequences={};sequence_counts=collections.Counter()
with zipfile.ZipFile(TABLE) as z:
    with z.open('tables/stable2.txt.gz') as f:
        with gzip.open(f,'rt') as g:
            for r in csv.DictReader(g,delimiter='\t'):
                if r['variant'] in wanted:sequences[r['variant']]=r;sequence_counts[r['variant']]+=1
lut=np.zeros(256,dtype='uint8');lut[[ord(x) for x in 'ACGT']]=[7,8,9,10]
def tok(s):return lut[np.frombuffer(s.encode(),dtype='uint8')]
excluded=collections.Counter();rows=[];ref=[];alt=[];activity=[];effect=[]
prior=collections.defaultdict(list)
with (ROOT/'outputs/stage5_score_queue.csv').open() as f:
    for r in csv.DictReader(f):prior[r['chr']].append(int(r['position_GRCh38']))
for c in prior:prior[c].sort()
for r in raw:
    key='chr'+r['SNP'];parts=key.split(':')
    if counts[r['SNP']]!=1:excluded['duplicate_assay_variant']+=1;continue
    if len(parts)!=4 or not parts[0].removeprefix('chr').isdigit() or not 1<=int(parts[0].removeprefix('chr'))<=22 or len(parts[2])!=1 or len(parts[3])!=1:excluded['not_autosomal_SNV']+=1;continue
    if sequence_counts[key]!=1 or key not in sequences:excluded['missing_or_duplicate_sequence_record']+=1;continue
    s=sequences[key]
    if s['duplicate']!='no':excluded['source_duplicate_flag']+=1;continue
    a,b=s['allele1_oligo'].upper(),s['allele2_oligo'].upper()
    if len(a)!=200 or len(b)!=200 or set(a+b)-set('ACGT'):excluded['unusable_construct']+=1;continue
    changes=[i for i,(x,y) in enumerate(zip(a,b)) if x!=y]
    if len(changes)!=1:excluded['not_single_substitution']+=1;continue
    p=changes[0]
    if p not in [99,100]:excluded['noncentered_SNV']+=1;continue
    if (a[p],b[p])!=(parts[2],parts[3]):excluded['source_allele_mismatch']+=1;continue
    try:
        aa,bb,ap,bp,skew=(float(r[k]) for k in ['A_log2FC','B_log2FC','A_logPadj_BF','B_logPadj_BF','Log2Skew'])
        assert all(math.isfinite(v) for v in [aa,bb,skew]) and all(not math.isnan(v) and v>=0 for v in [ap,bp])
    except (ValueError,AssertionError):excluded['unassessable_activity_or_effect']+=1;continue
    active=(abs(aa)>=1 and ap>=2) or (abs(bb)>=1 and bp>=2)
    if active:
        try:
            q=float(r['Skew_logFDR_act']);assert not math.isnan(q) and q>=0
        except (ValueError,AssertionError):excluded['active_but_unassessed_skew']+=1;continue
        label=int(q>=1)
    else:label=0
    chrom=parts[0];pos=int(parts[1]);mapped=out37.convert_coordinate(chrom,pos-1)
    if len(mapped)!=1 or mapped[0][0]!=chrom:excluded['unresolved_or_cross_chromosome_liftover']+=1;continue
    cc,pp,strand,score=mapped[0];reverse=back.convert_coordinate(cc,pp)
    if len(reverse)!=1 or reverse[0][:3]!=(chrom,pos-1,strand):excluded['nonreciprocal_mapping']+=1;continue
    ga,gb=(parts[2],parts[3]) if strand=='+' else (rc(parts[2]),rc(parts[3]));chrom_num=int(chrom[3:]);sp=split(chrom_num)
    if sp=='test':
        known=prior[cc];j=bisect.bisect_left(known,pp+1)
        if any(0<=k<len(known) and abs(known[k]-(pp+1))<=1000 for k in [j-1,j]):excluded['test_within_1kb_previously_inspected_project_queue']+=1;continue
    rows.append({'variant':key,'variant_GRCh38':f'{cc}:{pp+1}:{ga}:{gb}','study':'Siraj_OL42_SKNSH','chromosome':chrom_num,'position_hg38':pp+1,'genome_reference_allele':ga,'tested_alternate_allele':gb,'split':sp,'label':label,'active_screen':int(active),'changed_position':p,'source_row_id':r['ID']})
    ref.append(tok(a));alt.append(tok(b));activity.append((aa+bb)/2);effect.append(abs(skew))
print('Siraj paired/mapping eligibility passed',len(rows),'exclusions',dict(excluded),flush=True)
# Materialize Salomon inputs from their original source and v1 audited manifest.
with gzip.open(ROOT/'work/salomon_neuronal_variant_readout.tsv.gz','rt') as f:sal={r['chr_pos_ref_alt']:r for r in csv.DictReader(f,delimiter='\t')}
with np.load(OLD/'prepared.npz') as z:sdata={k:z[k] for k in z.files}
sm=list(csv.DictReader((ROOT/'outputs/sleep_frailty_ai_v1/split_manifest.csv').open()));sal_rows=[];sa=[];se=[]
for i,r in enumerate(sm):
    source=sal[r['variant']];c,pos,a,b=r['variant'].split('-');c=int(c)
    sal_rows.append({'variant':r['variant'],'variant_GRCh38':f'chr{c}:{pos}:{a}:{b}','study':'Salomon_NGN2','chromosome':c,'position_hg38':int(pos),'genome_reference_allele':a,'tested_alternate_allele':b,'split':split(c),'label':int(r['label']),'active_screen':int(source['NGN2_ref_significant_region']=='True' or source['NGN2_alt_significant_region']=='True'),'changed_position':int(r['changed_position']),'source_row_id':r['variant']})
    sa.append(.5*(float(source['NGN2_ref_region_logFC'])+float(source['NGN2_alt_region_logFC'])));se.append(abs(float(source['NGN2_variant_logFC'])))
alphabet=np.frombuffer(b'ACGT',dtype='uint8');byte_comp=bytes.maketrans(b'ACGT',b'TGCA')
def windows(x,p):
    # Salomon substitutions can be near fragment boundaries. Never invent an
    # empty centered window or treat short windows as equal 200-base contexts.
    dna=alphabet[x-7].tobytes();lo=max(0,p-199);hi=min(p,len(dna)-200)
    assert lo<=hi and 0<=p<len(dna)
    for start in range(lo,hi+1):
        core=dna[start:start+200];assert len(core)==200
        yield min(core,core.translate(byte_comp)[::-1])
all_rows=rows+sal_rows;all_ref=ref+list(sdata['ref']);all_alt=alt+list(sdata['alt']);sets={k:set() for k in ['validation','test']}
for r,a,b in zip(all_rows,all_ref,all_alt):
    if r['split'] in sets:
        p=r['changed_position'];sets[r['split']].update(windows(a,p));sets[r['split']].update(windows(b,p))
print('Verified mapped pairs and exact 200-base windows',len(rows),len(sal_rows),flush=True)
later_sets={'train':sets['test'] | sets['validation'],'validation':sets['test'],'test':set()}
keep=[]
for r,a,b in zip(all_rows,all_ref,all_alt):
    later=later_sets[r['split']]
    p=r['changed_position'];ok=r['split']=='test' or not (any(x in later for x in windows(a,p)) or any(x in later for x in windows(b,p)))
    keep.append(ok)
    if not ok:excluded[r['study']+'_exact_core_overlap_later_split']+=1
siraj_n=len(rows);siraj_keep=np.array(keep[:siraj_n]);sal_keep=np.array(keep[siraj_n:]);siraj_rows=[r for r,k in zip(rows,siraj_keep) if k];sal_rows=[r for r,k in zip(sal_rows,sal_keep) if k]
np.savez_compressed(WORK/'siraj_sequences_local.npz',ref=np.stack(ref)[siraj_keep],alt=np.stack(alt)[siraj_keep],position=np.array([r['changed_position'] for r in rows],dtype='int32')[siraj_keep],label=np.array([r['label'] for r in rows],dtype='float32')[siraj_keep],chromosome=np.array([r['chromosome'] for r in rows],dtype='int32')[siraj_keep],activity_mean=np.array(activity,dtype='float32')[siraj_keep],absolute_effect=np.array(effect,dtype='float32')[siraj_keep])
np.savez_compressed(WORK/'salomon_sequences_local.npz',ref=sdata['ref'][sal_keep],alt=sdata['alt'][sal_keep],position=sdata['position'][sal_keep],label=sdata['label'][sal_keep],chromosome=sdata['chromosome'][sal_keep],activity_mean=np.array(sa,dtype='float32')[sal_keep],absolute_effect=np.array(se,dtype='float32')[sal_keep],v1_row_index=np.flatnonzero(sal_keep))
for name,rr in [('siraj_manifest_local.csv',siraj_rows),('salomon_manifest_local.csv',sal_rows)]:
    with (WORK/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rr[0]));w.writeheader();w.writerows(rr)
combined=siraj_rows+sal_rows
with (WORK/'avi_query_pairs.csv').open('w',newline='') as f:
    fields=['variant','chromosome','position_hg38','genome_reference_allele','tested_alternate_allele'];w=csv.DictWriter(f,fieldnames=fields);w.writeheader();seen=set()
    for r in combined:
        if r['variant_GRCh38'] in seen:continue
        seen.add(r['variant_GRCh38']);w.writerow({'variant':r['variant_GRCh38'],'chromosome':'chr'+str(r['chromosome']),**{k:r[k] for k in fields[2:]}})
summary={'source_files':{'Siraj_fixed_library':{'path':str(SOURCE),'sha256':sha(SOURCE)},'Siraj_sequences_archive':{'path':str(TABLE),'sha256':sha(TABLE)},'Salomon_readout_sha256':sha(ROOT/'work/salomon_neuronal_variant_readout.tsv.gz')},'source_license':'Siraj Zenodo record cc-by-nc-4.0; research use and attribution; do not imply clinical use','raw_fixed_library_center_rows':len(raw),'excluded':dict(excluded),'studies':{},'AVI_unique_allele_queries':len(seen),'split_rule':'Chromosomes 1-3 test, 4-5 validation, 6-22 train across both studies','test_status':'No v3 model predictions evaluated. Salomon study previously developed; Siraj held-chromosome model results not yet inspected. Previously inspected project queue neighborhoods excluded from Siraj test.','auxiliary_targets':'Measured activity mean and absolute allelic effect are training targets only; never inference inputs','canonical_core_leakage_check':'Exact canonical 200-base reference/alternate windows containing the tested substitution, including all valid windows in 270-base Salomon fragments. Remove earlier-split matches to later split. No empty/short windows.','limitations':['Fixed source-row screen, not authors final best-library hit list','Different neuron states, construct lengths, statistical detection and study selection','Source sequence and original genome pretraining/model benchmark exposure not independently excluded','No participant-level sleep/frailty prediction or genetic-score validation']}
for study,rr in [('Siraj_OL42_SKNSH',siraj_rows),('Salomon_NGN2',sal_rows)]:summary['studies'][study]={sp:{'n':sum(r['split']==sp for r in rr),'detected_calls':sum(r['label'] for r in rr if r['split']==sp)} for sp in ['train','validation','test']}
(OUT/'data_audit.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2),flush=True)
