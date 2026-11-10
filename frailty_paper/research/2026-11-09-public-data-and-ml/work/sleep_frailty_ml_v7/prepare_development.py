"""Magnitude targets with separate units, complete-study exclusion and shared cell loci."""
import gzip,collections,bisect
from common import *
from reference import get_reference
from prepare_reserved_design import json_list,rc
GROUP={'neuronal_induced':'induced_excitatory_neuron','neuronal_tumor':'neuronal_tumor','kidney':'kidney_derived','liver':'hepatocyte_tumor','colon':'colon_tumor','blood':'hematopoietic','glioma':'glioma','cervical':'cervical','skin':'skin','melanoma':'melanoma','prostate':'prostate'}
def table(path):
    op=gzip.open if str(path).endswith(('.gz',)) else open
    with op(path,'rt',newline='') as f:return list(csv.DictReader(f,delimiter='\t' if str(path).endswith('.gz') and 'sequence_bank' not in str(path) and 'original_effects' not in str(path) and 'continuous_effects' not in str(path) else ','))
def qc_counts(path):
    out=collections.defaultdict(dict)
    with gzip.open(path,'rt') as f:
        for r in csv.DictReader(f,delimiter='\t'):
            assert r['replicate'] not in out[r['oligo_name']]
            out[r['oligo_name']][r['replicate']]=[float(r[k]) if r.get(k) not in ['',None,'NA','NaN'] else np.nan for k in ['n_bc','dna_counts','rna_counts']]
    return out
def source_effects(path):
    with gzip.open(path,'rt') as f:rows=list(csv.DictReader(f,delimiter='\t'))
    counts=collections.Counter(r['variant_id'] for r in rows)
    return {r['variant_id']:r for r in rows if counts[r['variant_id']]==1}
def valid_qc(qc,ref,alt,minimum=15):
    a,b=qc.get(ref),qc.get(alt)
    if a is None or b is None or set(a)!=set(b) or len(a)<3:return False
    v=np.array([x for r in [a,b] for x in r.values()]);return bool(np.isfinite(v).all() and np.all(v[:,0]>=minimum) and np.all(v[:,1]>0) and np.all(v[:,2]>=0))
def main():
    assert not (OUT/'model_search_started.json').exists()
    prior=assert_previous_intact();audit=json.loads((OUT/'reserved_design_and_original_source_audit.json').read_text());assert sha(DATA/'reserved_physical_design.csv')==audit['physical_design_sha256']
    test=list(csv.DictReader((DATA/'reserved_physical_design.csv').open()));points=collections.defaultdict(set);test_cores=set()
    for r in test:
        points[r['chromosome']].add(int(r['position_hg38']))
        for s in [r['reference'],r['alternate']]:
            p=int(r['changed_position']);lo=max(0,p-149);hi=min(p,len(s)-150)
            for j in range(lo,hi+1):v=s[j:j+150];test_cores.add(v);test_cores.add(rc(v))
    points={c:sorted(v) for c,v in points.items()};prefixes={s[:16] for s in test_cores}
    def matches_core(r):
        p=int(r['changed_position']);a=r['reference'];lo=max(0,p-149);hi=min(p,len(a)-150)
        for s in [a,r['alternate']]:
            for j in range(lo,hi+1):
                if s[j:j+16] in prefixes and s[j:j+150] in test_cores:return True
        return False
    with gzip.open(V6/'sequence_bank.csv.gz','rt') as f:oldbank=list(csv.DictReader(f))
    bank=[];obs=[];source_exclusions=collections.Counter()
    def add_bank(r,feature_source,feature_index):
        i=len(bank);bank.append({'bank_id':i,'variant_GRCh38':r['variant_GRCh38'],'chromosome':r['chromosome'],'position_hg38':int(r['position_hg38']),'reference':r['reference'],'alternate':r['alternate'],'changed_position':int(r['changed_position']),'context_feature_source':feature_source,'context_feature_index':int(feature_index)});return i
    def observation(i,study,experiment,cell,group,value,unit='absolute_log2_effect',classification='',parent=''):
        value=float(value)
        if not np.isfinite(value) or value<0:source_exclusions[study+'_unassessed_magnitude']+=1;return
        obs.append({'bank_id':i,'study':study,'experiment':experiment,'cell_line':cell,'cell_group':group,'absolute_effect':value,'effect_unit':unit,'source_detection':classification,'parent_locus':parent})
    for r in oldbank:assert add_bank(r,'v6',r['bank_id'])==int(r['bank_id'])
    with (V6/'measured_observations.csv').open() as f:
        for r in csv.DictReader(f):observation(int(r['bank_id']),r['study'],r['experiment'],r['cell'],GROUP[r['cell_group']],r['absolute_effect'],'saturation_coefficient' if r['study']=='Kircher2019' else 'absolute_log2_effect',r['label'],r.get('parent_locus',''))
    print('Existing measured data materialized',len(bank),len(obs),flush=True)
    # Recovered neuronal outcomes are explicitly observed V7 development data.
    from_old=json.loads((V6OUT/'kosicki_continuous_recovery_audit.json').read_text());p=V6OUT/'kosicki_reconstructed_continuous_effects.csv';assert sha(p)==from_old['reconstructed_effects_sha256']
    sf=V6/'sources/kosicki_ENCFF752NFF.fasta.gz';sequences={};name=None;pieces=[]
    with gzip.open(sf,'rt') as f:
        for line in f:
            if line.startswith('>'):
                if name is not None:sequences[name]=''.join(pieces).upper()
                name=line[1:].strip();pieces=[]
            else:pieces.append(line.strip())
    if name is not None:sequences[name]=''.join(pieces).upper()
    with p.open() as f:
        for r in csv.DictReader(f):
            r={**r,'reference':sequences[r['reference_name']],'alternate':sequences[r['name']]};i=add_bank(r,'kosicki',r['design_index']);observation(i,'Kosicki2025','NGN2_lentiMPRA','WTC11_NGN2','induced_excitatory_neuron',abs(float(r['mean_log2_ratio_difference'])),'absolute_log2_effect','',r['reference_name'])
    del sequences,oldbank
    # The prior held-out HepG2 library is now development. Its exact design
    # and identical SNPs in HEK293T share one sequence row and one locus role.
    hepatic=source_effects(V6/'sources/recovery_IGVFFI4642CBGA.gz');kidney=source_effects(SOURCES/'IGVFFI0588HACZ.gz');kidney_qc=qc_counts(SOURCES/'IGVFFI9432KCNV.gz')
    with (V6/'additional_study_blind_design.csv').open() as f:paired=list(csv.DictReader(f))
    for j,r in enumerate(paired):
        i=add_bank(r,'hepg2',j);s=hepatic.get(r['SPDI'])
        if s is not None:observation(i,'IGVF240K','240K_HepG2','HepG2','hepatocyte_tumor',abs(float(s['log2FoldChange'])))
        s=kidney.get(r['SPDI'])
        if s is None:source_exclusions['HEK293T_no_unique_reported_effect']+=1;continue
        if not valid_qc(kidney_qc,r['reference_name'],r['name']):source_exclusions['HEK293T_pair_failed3rep15barcode_QC']+=1;continue
        observation(i,'IGVF240K','240K_HEK293T','HEK293T','kidney_derived',abs(float(s['log2FoldChange'])))
    del hepatic,kidney,kidney_qc,paired
    print('Observed neuronal and matched kidney/liver development added',len(bank),len(obs),flush=True)
    # New NPC source. Source-paper culture, not the inconsistent portal
    # targeted-sample label, defines the explicit cell context.
    design_path=SOURCES/'IGVFFI4914OUJH.tsv.gz'
    with gzip.open(design_path,'rt') as f:raw=list(csv.DictReader(f,delimiter='\t'))
    refs=collections.defaultdict(list);alts=[]
    for r in raw:
        if r['category']!='variant' or r['class']!='test' or r['ref']!='GRCh38':continue
        try:sp,aa,vc,ps=map(json_list,[r['SPDI'],r['allele'],r['variant_class'],r['variant_pos']])
        except (ValueError,TypeError,AssertionError):source_exclusions['ASD_missing_structured_variant_annotation']+=1;continue
        if len(sp)!=1 or vc!=['SNV'] or len(ps)!=1 or aa not in [['ref'],['alt']]:source_exclusions['ASD_not_single_SNV']+=1;continue
        r={**r,'SPDI':sp[0],'changed_position':int(ps[0])}
        if aa==['ref']:refs[r['SPDI']].append(r)
        else:alts.append(r)
    genome=get_reference();effects=source_effects(SOURCES/'IGVFFI7498KQLC.gz');qcs=qc_counts(SOURCES/'IGVFFI3167QVNF.gz');numbers=collections.Counter(r['SPDI'] for r in alts)
    for b in alts:
        spdi=b['SPDI'];rr=refs.get(spdi,[])
        if len(rr)!=1 or numbers[spdi]!=1:source_exclusions['ASD_nonunique_physical_pair']+=1;continue
        a=rr[0];ref,alt=a['sequence'].upper(),b['sequence'].upper();ps=[k for k,(x,y) in enumerate(zip(ref,alt)) if x!=y]
        if len(ref)!=len(alt) or len(ref)!=200 or set(ref+alt)-set('ACGT') or len(ps)!=1:source_exclusions['ASD_invalid200bp_physical_SNV']+=1;continue
        j=ps[0];acc,pp,ra,aa=spdi.split(':');pp=int(pp);start,end=int(a['start']),int(a['end']);strand=a['strand']
        expected=start+j if strand=='+' else end-j-1;alleles=(ref[j],alt[j]) if strand=='+' else (rc(ref[j]),rc(alt[j]))
        if any(a[k]!=b[k] for k in ['chr','start','end','strand']) or j!=a['changed_position'] or j!=b['changed_position'] or end-start!=200 or expected!=pp or alleles!=(ra,aa):source_exclusions['ASD_annotation_or_allele_disagreement']+=1;continue
        if genome.sequence(a['chr'],start,end)!=(ref if strand=='+' else rc(ref)):source_exclusions['ASD_full_reference_genome_mismatch']+=1;continue
        an,bn=a['name'].replace(' ','_'),b['name'].replace(' ','_')
        if not valid_qc(qcs,an,bn):source_exclusions['ASD_pair_failed3rep15barcode_QC']+=1;continue
        s=effects.get(spdi)
        if s is None:source_exclusions['ASD_no_unique_reported_effect']+=1;continue
        r={'variant_GRCh38':f"{a['chr']}:{pp+1}:{ra}:{aa}",'chromosome':a['chr'],'position_hg38':pp+1,'reference':ref,'alternate':alt,'changed_position':j};i=add_bank(r,'new',-1);observation(i,'Koesterich2023','ASD_N2_lentiMPRA','H1_N2_NPC','neural_progenitor',abs(float(s['log2FoldChange'])))
    del raw,refs,alts,effects,qcs
    print('NPC physical sequence and measured magnitude audit complete',len(bank),len(obs),flush=True)
    # Connect before exclusions so deleting a test-overlapping bridge cannot
    # silently split its neighboring training component into allowed parts.
    parent=list(range(len(bank)))
    def find(i):
        while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
        return i
    def union(a,b):
        a,b=find(a),find(b)
        if a!=b:parent[max(a,b)]=min(a,b)
    chrom=collections.defaultdict(list);core_seen={}
    for r in bank:
        i=r['bank_id'];chrom[r['chromosome']].append(i);p=int(r['changed_position']);a=r['reference'];start=max(0,min(p-75,len(a)-150));s=a[start:start+150];s=min(s,rc(s))
        if s in core_seen:union(i,core_seen[s])
        else:core_seen[s]=i
    for ids in chrom.values():
        ids.sort(key=lambda i:bank[i]['position_hg38'])
        for a,b in zip(ids,ids[1:]):
            if bank[b]['position_hg38']-bank[a]['position_hg38']<=2000:union(a,b)
    bad_roots=set();overlap=collections.Counter()
    for i,r in enumerate(bank):
        ps=points.get(r['chromosome'],[]);p=r['position_hg38'];j=bisect.bisect_left(ps,p)
        if any(0<=k<len(ps) and abs(ps[k]-p)<=2000 for k in [j-1,j]):bad_roots.add(find(i));overlap['within2kb_reserved_position']+=1
        elif matches_core(r):bad_roots.add(find(i));overlap['exact_reserved_variant_containing150bp_core']+=1
        if i%50000==0:print('Reserved whole-locus/core exclusion audit',i,'/',len(bank),flush=True)
    roots=collections.defaultdict(list)
    for r in bank:roots[find(r['bank_id'])].append(r['bank_id'])
    allowed={i for root,ids in roots.items() if root not in bad_roots for i in ids};old_to_new={old:new for new,old in enumerate(sorted(allowed))};newbank=[]
    for old in sorted(allowed):newbank.append({**bank[old],'bank_id':old_to_new[old]})
    for root,ids in roots.items():
        if root in bad_roots:continue
        coords=sorted({(bank[i]['chromosome'],bank[i]['position_hg38']) for i in ids});h=hashlib.sha256(('V7_effect_cell_locus_20261107|'+'|'.join(f'{c}:{p}' for c,p in coords)).encode()).hexdigest();fold=int(h[:8],16)%10;role='calibration' if fold<2 else 'validation' if fold<4 else 'train'
        for i in ids:newbank[old_to_new[i]].update({'locus_group':h[:24],'split':role})
    accepted=[];excluded_observations=collections.Counter()
    for r in obs:
        old=r['bank_id']
        if old not in allowed:excluded_observations[r['study']+'|'+r['cell_group']]+=1;continue
        b=newbank[old_to_new[old]];accepted.append({**r,'bank_id':b['bank_id'],'locus_group':b['locus_group'],'split':b['split']})
    del bank,obs,core_seen,parent
    with gzip.open(DATA/'sequence_bank.csv.gz','wt',newline='') as f:w=csv.DictWriter(f,fieldnames=list(newbank[0]));w.writeheader();w.writerows(newbank)
    with (DATA/'measured_magnitudes.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(accepted[0]));w.writeheader();w.writerows(accepted)
    inventory=collections.defaultdict(collections.Counter)
    for r in accepted:inventory[r['study']+'|'+r['cell_group']+'|'+r['effect_unit']][r['split']]+=1
    report={'observations':len(accepted),'paired_sequences':len(newbank),'primary_log2_observations':sum(r['effect_unit']=='absolute_log2_effect' for r in accepted),'auxiliary_saturation_coefficient_observations':sum(r['effect_unit']=='saturation_coefficient' for r in accepted),'inventory':{k:dict(v) for k,v in inventory.items()},'source_eligibility_exclusions':dict(source_exclusions),'test_overlap_triggers':dict(overlap),'excluded_whole_locus_components':len(bad_roots),'excluded_observations_by_source_cell':dict(excluded_observations),'partition':'All source/cell repeats together; <=2kb connected components and shared central150bp reference contexts. Any component touching reserved positions within2kb or an exact reserved150bpvariant-containing core is entirely excluded. SHA256 seed V7_effect_cell_locus_20261107;20%calibration,20%selection,60%training.','reserved_effect_values_used_for_training_selection_or_calibration':False,'observed_prior_test_role':'Kosicki2025 and IGVF240K HepG2 are already observed V7 development; never fresh tests again.','units':'Primary absolute log2reporter effect. Kircher coefficients explicitly separate auxiliary units, never silently treated as log2expression.','NPC_provenance':'Koesterich2023 section4.7.1 describes H1-derived N2 early neural progenitors. Portal glutamatergic-neuron7day label is inconsistent and not used as a verified cell-state claim.','V6_lock_sha256':prior,'source_assignment_plan_sha256':sha(WORK/'SOURCE_ASSIGNMENT_AMENDMENT.txt'),'reserved_design_sha256':sha(DATA/'reserved_physical_design.csv'),'outputs':{n:sha(DATA/n) for n in ['sequence_bank.csv.gz','measured_magnitudes.csv']}}
    write_json(OUT/'development_data_audit.json',report);print('V7 MAGNITUDE DATA AND WHOLE-STUDY EXCLUSIONS READY',len(accepted),'observations;',len(newbank),'paired DNA sequences',flush=True)
if __name__=='__main__':main()
