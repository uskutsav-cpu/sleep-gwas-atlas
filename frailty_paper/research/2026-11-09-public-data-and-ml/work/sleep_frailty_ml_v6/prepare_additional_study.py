"""Audit an additional study's sequence/QC design without model-outcome comparisons."""
import gzip,collections,bisect,re
from common import *
from additional_reference import get_reference
from prepare_blind_test_design import rc
from prepare_measured_data import core_hashes
def main():
    assert not (OUT/'additional_test_evaluation_started.json').exists()
    original=json.loads((OUT/'final_model_lock.json').read_text())
    for path,h in original['files'].items():assert sha(path)==h,path
    excluded=collections.Counter();qc=collections.defaultdict(dict)
    with gzip.open(SOURCES/'recovery_IGVFFI1642PFRO.gz','rt') as f:
        for r in csv.DictReader(f,delimiter='\t'):
            rep=int(r['replicate']);assert rep in [1,2,3] and rep not in qc[r['oligo_name']]
            # Read counts for assessability only, not effect labels/rank comparisons.
            qc[r['oligo_name']][rep]=[float(r[k]) for k in ['n_bc','dna_counts','rna_counts']]
    design=[]
    with gzip.open(SOURCES/'recovery_IGVFFI2260ZLFH.tsv.gz','rt') as f:
        for r in csv.DictReader(f,delimiter='\t'):
            if r['category']!='variant' or r['class']!='test' or r['ref']!='GRCh38' or r['name'] not in qc:continue
            r={**r,'spdis':json.loads(r['SPDI']),'alleles':json.loads(r['allele']),'positions':json.loads(r['variant_pos'])};design.append(r)
    refmap=collections.defaultdict(list)
    for r in design:
        if all(a=='ref' for a in r['alleles']):
            for spdi in r['spdis']:refmap[(spdi,r['chr'],r['start'],r['end'],r['strand'])].append(r)
    candidates=collections.defaultdict(list)
    for r in design:
        if r['alleles']!=['alt']:continue
        if len(r['spdis'])!=1 or json.loads(r['variant_class'])!=['SNV']:excluded['not_single_SNV']+=1;continue
        spdi=r['spdis'][0];refs=refmap.get((spdi,r['chr'],r['start'],r['end'],r['strand']),[])
        if len(refs)!=1:excluded['missing_or_multiple_paired_reference_designs']+=1;continue
        a=refs[0];ref,alt=a['sequence'].upper(),r['sequence'].upper();changes=[i for i,(x,y) in enumerate(zip(ref,alt)) if x!=y]
        if len(ref)!=200 or len(alt)!=200 or set(ref+alt)-set('ACGT') or len(changes)!=1:excluded['invalid200bp_single_substitution']+=1;continue
        pos=changes[0];accession,p0,ra,aa=spdi.split(':');p0=int(p0)
        if len(ra)!=1 or len(aa)!=1 or r['positions']!=[pos] or r['chr']!='chr'+{23:'X',24:'Y'}.get(int(accession.split('.')[0].split('_')[1]),str(int(accession.split('.')[0].split('_')[1]))):excluded['SPDI_metadata_disagrees']+=1;continue
        start,end=int(r['start']),int(r['end']);expected=start+pos if r['strand']=='+' else end-pos-1
        rbase,abase=(ref[pos],alt[pos]) if r['strand']=='+' else (rc(ref[pos]),rc(alt[pos]))
        if expected!=p0 or (rbase,abase)!=(ra,aa) or end-start!=200:excluded['coordinate_or_allele_orientation_disagrees']+=1;continue
        candidates[spdi].append({'name':r['name'],'reference_name':a['name'],'SPDI':spdi,'variant_GRCh38':f"{r['chr']}:{p0+1}:{ra}:{aa}",'chromosome':r['chr'],'position_hg38':p0+1,'start':start,'end':end,'strand':r['strand'],'reference':ref,'alternate':alt,'changed_position':pos,'source_class':r['class']})
    pairs=[]
    for spdi,items in candidates.items():
        if len(items)!=1:excluded['multiple_constructs_for_one_source_aggregated_variant']+=len(items);continue
        r=items[0];values=[qc[n] for n in [r['name'],r['reference_name']]]
        if any(set(v)!={1,2,3} for v in values):excluded['not_all3technical_replicates']+=1;continue
        counts=np.array([z for v in values for z in v.values()])
        if not np.isfinite(counts).all() or np.any(counts[:,0]<15) or np.any(counts[:,1]<=0) or np.any(counts[:,2]<0):excluded['failed_prespecified_barcode_or_count_QC']+=1;continue
        r['minimum_barcode_count']=float(counts[:,0].min());pairs.append(r)
    # Remove every sequence used anywhere in fitting OR confidence calibration.
    with gzip.open(DATA/'sequence_bank.csv.gz','rt') as f:bank=list(csv.DictReader(f))
    points=collections.defaultdict(list);alleles=set()
    for r in bank:points[r['chromosome']].append(int(r['position_hg38']));alleles.add(r['variant_GRCh38'])
    points={c:sorted(set(p)) for c,p in points.items()};far=[]
    for r in pairs:
        p=int(r['position_hg38']);ps=points.get(r['chromosome'],[]);j=bisect.bisect_left(ps,p)
        if r['variant_GRCh38'] in alleles:excluded['exact_previously_used_allele']+=1;continue
        if any(0<=k<len(ps) and abs(ps[k]-p)<=1000 for k in [j-1,j]):excluded['within1kb_previously_used_variant']+=1;continue
        far.append(r)
    cores=set()
    for i,r in enumerate(bank):
        cores.update(core_hashes(r['reference'],r['alternate'],int(r['changed_position'])))
        if i%30000==0:print('Additional study overlap dictionary',i,'/',len(bank),flush=True)
    disjoint=[]
    for r in far:
        if any(h in cores for h in core_hashes(r['reference'],r['alternate'],int(r['changed_position']))):excluded['exact150bp_variant_containing_core']+=1;continue
        disjoint.append(r)
    del cores,bank
    genome=get_reference();kept=[]
    for r in disjoint:
        expected=genome.sequence(r['chromosome'],r['start'],r['end']);actual=r['reference'] if r['strand']=='+' else rc(r['reference'])
        if expected!=actual:excluded['full_construct_reference_genome_mismatch']+=1;continue
        kept.append(r)
    assert len(kept)>1000 and len({r['SPDI'] for r in kept})==len(kept)
    p=DATA/'additional_study_blind_design.csv'
    with p.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(kept[0]));w.writeheader();w.writerows(kept)
    report={'source':'IGVFDS6149YGOL HepG2 assay / IGVFDS1589ODOW library','parent_sample':'IGVFSM9009DVDG','source_build':'GRCh38,SPDI0based positions converted to1based genomic variants','entire_assay_oligo_names':len(qc),'test_SNV_variant_ids_before_unique_construct_QC':len(candidates),'single_construct_QC_passed_pairs':len(pairs),'eligible_strict_nonoverlapping_pairs':len(kept),'excluded':dict(excluded),'entire_reference_construct_verified':True,'source_original_model_lock_sha256':sha(OUT/'final_model_lock.json'),'design_sha256':sha(p),'source_files':{n:sha(SOURCES/n) for n in ['recovery_IGVFFI1642PFRO.gz','recovery_IGVFFI2260ZLFH.tsv.gz','recovery_IGVFDS6149YGOL.json']},'model_effect_comparison_performed':False,'scope':'New measurement study, same established HepG2 cell line;3technical replicates, not new independent donors. Source variants were CNN-selected, not a population-random variant sample. Source processed effects remain sealed from prediction comparisons.'}
    write_json(OUT/'additional_study_blind_design_audit.json',report);print('ADDITIONAL STUDY DESIGN CHECKED',len(kept),'pairs',dict(excluded),flush=True)
if __name__=='__main__':main()
