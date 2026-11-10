"""Conservatively exclude every original saturation locus, covering CAGI5."""
import bisect,csv,hashlib,json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
DEST=ROOT/'outputs/sleep_frailty_ml_bounded/AVI_exposure'
PREVIOUS=ROOT/'outputs/sleep_frailty_ml_native_stable/benchmark_exposure'


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    source=DEST/'Kircher_all_author_saturation_variants.tsv'
    receipt=json.loads(source.with_suffix('.retrieval.json').read_text());assert sha(source)==receipt['sha256']
    positions={};elements=set();nrows=0
    for r in csv.DictReader(source.open(),delimiter='\t'):
        positions.setdefault('chr'+r['Chromosome'].removeprefix('chr'),set()).add(int(r['Position']))
        elements.add(r['Element']);nrows+=1
    positions={c:sorted(p) for c,p in positions.items()}
    # The 2025 PDF figure was visually checked; the 2026 published HTML
    # independently confirms CAGI5 saturation MPRA as an evaluated benchmark.
    html=DEST/'AlphaGenome_base_2026_published.html'
    assert b'CAGI5' in html.read_bytes() and b'Saturation mutagenesis of twenty' in html.read_bytes()
    profile=ROOT/'outputs/sleep_frailty_ml_channels/source_candidates/Ulirsch_complete_exposure_profile.csv'
    allrows=list(csv.DictReader(profile.open()));excluded=set();hits=[]
    for r in allrows:
        chrom,pos,*_=r['variant_GRCh38'].split(':');pos=int(pos);values=positions.get(chrom,[])
        j=bisect.bisect_left(values,pos);nearest=values[max(j-1,0):min(j+1,len(values))]
        if nearest and min(abs(pos-p) for p in nearest)<=2145:
            excluded.add(r['joint_component_after_all_source_DNA'])
            hits.append({'variant':r['variant_GRCh38'],'component':r['joint_component_after_all_source_DNA']})
    pools=[('base',PREVIOUS/'Ulirsch_candidates_after_all_documented_reader_benchmarks.csv'),
           ('matched_AVI',DEST/'Ulirsch_candidates_outside_AVI_development_chromosomes.csv')]
    results={}
    for name,pool in pools:
        rows=list(csv.DictReader(pool.open()));kept=[r for r in rows if r['joint_component_after_all_source_DNA'] not in excluded]
        out=DEST/f'Ulirsch_{name}_candidates_after_CAGI5_guard.csv'
        with out.open('w') as handle:
            writer=csv.DictWriter(handle,fieldnames=list(allrows[0]));writer.writeheader();writer.writerows(kept)
        results[name]={'input_sha256':sha(pool),'candidate_tests_before':len(rows),'candidate_tests_after':len(kept),
           'candidate_variants_after':len({r['variant_GRCh38'] for r in kept}),
           'candidate_components_after':len({r['joint_component_after_all_source_DNA'] for r in kept}),
           'output':str(out.relative_to(ROOT)),'output_sha256':sha(out)}
    result={'passed':True,'original_saturation_source_sha256':sha(source),'original_source_rows':nrows,
       'original_source_element_names':sorted(elements),'published_2026_paper_URL':'https://www.nature.com/articles/s41586-025-10014-0',
       'published_2026_HTML_sha256':sha(html),'2025_PDF_visual_inventory_checked':'Extended Data Figure 8, PDF page 77',
       'scope':'All released saturation element positions, including extra loci and indels beyond the published CAGI5 SNV subset, with 2145bp neighborhood and whole-source-component propagation.',
       'source_profile_sha256':sha(profile),'affected_full_source_components':len(excluded),'all_full_source_matches':hits,'candidate_pools':results,
       'new_study_outcome_files_read':False,'AVI_scores_retrieved':False,'model_performance_calculated':False,
       'limits':'Covers the documented saturation MPRA benchmark. Does not prove no prior genomic input, QTL knowledge, or other base-model exposure. These remain candidate identities before sealed count QC and model/confidence freezing.',
       'code_sha256':sha(Path(__file__))}
    (DEST/'base_AlphaGenome_CAGI5_exposure_receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(results),flush=True)


if __name__=='__main__':main()
