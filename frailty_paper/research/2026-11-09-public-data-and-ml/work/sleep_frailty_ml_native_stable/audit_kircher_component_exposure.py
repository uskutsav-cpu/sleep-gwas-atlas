"""Propagate original four-locus saturation benchmark exposure by component."""
import csv,hashlib,json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
DEST=ROOT/'outputs/sleep_frailty_ml_native_stable/benchmark_exposure'


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    retrieval=json.loads((DEST/'kircher/original_saturation_inventory_retrieval.json').read_text())
    loci=[]
    for record in retrieval['files']:
        path=DEST/'kircher'/record['name'];assert sha(path)==record['sha256']
        data=list(csv.DictReader(path.open(),delimiter='\t'))
        chrom={r['Chromosome'] for r in data};assert len(chrom)==1
        positions=[int(r['Position']) for r in data]
        loci.append({'file':record['name'],'chromosome':'chr'+next(iter(chrom)).removeprefix('chr'),
             'lowest_original_position':min(positions),'highest_original_position':max(positions),
             'rows':len(data),'sha256':record['sha256']})
    source=ROOT/'outputs/sleep_frailty_ml_channels/source_candidates/Ulirsch_complete_exposure_profile.csv'
    allrows=list(csv.DictReader(source.open()));touched=set();hits=[]
    for r in allrows:
        chrom,pos,*_=r['variant_GRCh38'].split(':');pos=int(pos)
        for locus in loci:
            # Both possible base-index conventions covered conservatively by
            # one extra bp. Whole original locus, not only significant variants.
            if chrom==locus['chromosome'] and locus['lowest_original_position']-2145<=pos<=locus['highest_original_position']+2145:
                component=r['joint_component_after_all_source_DNA'];touched.add(component)
                hits.append({'variant':r['variant_GRCh38'],'component':component,'benchmark':locus['file']})
    prior=DEST/'Ulirsch_candidates_after_author_benchmark_guard.csv';candidates=list(csv.DictReader(prior.open()))
    kept=[r for r in candidates if r['joint_component_after_all_source_DNA'] not in touched]
    assert not {r['joint_component_after_all_source_DNA'] for r in kept}&touched
    final=DEST/'Ulirsch_candidates_after_all_documented_reader_benchmarks.csv'
    with final.open('w') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(allrows[0]));writer.writeheader();writer.writerows(kept)
    receipt={'passed':True,'source_profile_sha256':sha(source),'previous_candidate_sha256':sha(prior),
       'published_saturation_benchmark_loci':loci,'all_matching_source_rows':hits,
       'whole_source_components_near_saturation_loci':len(touched),'candidate_designs_before':len(candidates),
       'candidate_designs_after':len(kept),'candidate_variants_after':len({r['variant_GRCh38'] for r in kept}),
       'candidate_components_after':len({r['joint_component_after_all_source_DNA'] for r in kept}),
       'scope':'All eight released experiments/replicates/orientations across the four published saturation benchmark loci; interval widened 2145bp to cover full fragment, neighborhood and any one-base convention ambiguity.',
       'new_study_outcome_files_read':False,'performance_calculated':False,
       'outside_validation_ready':False,'remaining':'Supplementary-methods exposure review; scientific candidate and confidence rules must be fixed before evaluating new-study outcomes.',
       'candidate_file_sha256':sha(final),'code_sha256':sha(Path(__file__))}
    (DEST/'saturation_benchmark_component_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:v for k,v in receipt.items() if k.startswith('candidate_') or k=='whole_source_components_near_saturation_loci'}))


if __name__=='__main__':main()
