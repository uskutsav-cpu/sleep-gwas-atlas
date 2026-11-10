"""Reserve an identity-only matched test pool outside AVI development chromosomes."""
import csv,hashlib,json,ssl,urllib.request,re,html
from pathlib import Path
import certifi

ROOT=Path(__file__).resolve().parents[2]
DEST=ROOT/'outputs/sleep_frailty_ml_bounded/AVI_exposure'
PREVIOUS=ROOT/'outputs/sleep_frailty_ml_native_stable/benchmark_exposure'
TRAIN=['1','4','7','8','10','13','15']
VALIDATION=['2','5','11','14','17','20','22','X']


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    DEST.mkdir(exist_ok=True)
    url='https://pmc.ncbi.nlm.nih.gov/articles/PMC13596691/'
    with urllib.request.urlopen(url,context=ssl.create_default_context(cafile=certifi.where()),timeout=30) as response:raw=response.read(2000001)
    assert len(raw)<=2000000
    text=html.unescape(re.sub('<[^>]+>',' ',raw.decode()));text=re.sub(r'\s+',' ',text)
    assert 'Data Split' in text and 'Training: proxy neutral' in text
    assert '1, 4, 7, 8, 10, 13, 15' in text and '2, 5, 11, 14, 17, 20, 22, X' in text
    (DEST/'Atlas_original_publication.html').write_bytes(raw)
    profile=ROOT/'outputs/sleep_frailty_ml_channels/source_candidates/Ulirsch_complete_exposure_profile.csv'
    allrows=list(csv.DictReader(profile.open()));bad={'chr'+c for c in TRAIN+VALIDATION}
    excluded_components={r['joint_component_after_all_source_DNA'] for r in allrows if r['chromosome'] in bad}
    prior=PREVIOUS/'Ulirsch_candidates_after_all_documented_reader_benchmarks.csv';rows=list(csv.DictReader(prior.open()))
    kept=[r for r in rows if r['joint_component_after_all_source_DNA'] not in excluded_components]
    assert all(r['chromosome'] not in bad for r in kept)
    target=DEST/'Ulirsch_candidates_outside_AVI_development_chromosomes.csv'
    with target.open('w') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(allrows[0]));writer.writeheader();writer.writerows(kept)
    receipt={'passed_publication_and_identity_screen':True,'publication':'https://pmc.ncbi.nlm.nih.gov/articles/PMC13596691/',
       'publication_type':'2026-09-20 version-1 preprint; not peer reviewed','public_primary_source_url':url,
       'original_HTML_sha256':hashlib.sha256(raw).hexdigest(),'AVI_proxy_label_training_chromosomes':TRAIN,
       'AVI_selection_validation_chromosomes':VALIDATION,'source_profile_sha256':sha(profile),
       'previous_candidate_sha256':sha(prior),'candidate_tests_before':len(rows),'candidate_tests_after':len(kept),
       'candidate_variants_after':len({r['variant_GRCh38'] for r in kept}),
       'candidate_components_after':len({r['joint_component_after_all_source_DNA'] for r in kept}),
       'eligibility_scope':'Exclude whole linked source components if any profile member lies on an AVI training or validation chromosome, conservatively excluding whole validation chromosomes rather than only known validation variants. Retain existing DNA/history/pretraining/author-benchmark exclusions.',
       'does_not_establish':'Does not prove absence from every base AlphaGenome genomic training or benchmark resource. A gnomAD rarity proxy is not an experimental effect-size target.',
       'future_comparison':'Compare both permitted-AVI and zero-AVI controls on exactly the same eligible cases, after model/confidence rules are fixed. This is a candidate identity pool before applying the already sealed count-QC rule; no outcome selection or accuracy calculation performed.',
       'new_study_outcomes_read':False,'AVI_scores_retrieved':False,'model_performance_calculated':False,
       'candidate_file_sha256':sha(target),'code_sha256':sha(Path(__file__))}
    (DEST/'AVI_development_chromosome_exposure_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:v for k,v in receipt.items() if k.startswith('candidate_')}),flush=True)


if __name__=='__main__':main()
