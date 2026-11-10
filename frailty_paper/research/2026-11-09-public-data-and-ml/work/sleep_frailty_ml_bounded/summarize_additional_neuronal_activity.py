"""Summarize additional physical constructs, parsing decimal-form coordinates."""
import collections,csv,gzip,json,math
from decimal import Decimal,InvalidOperation
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'outputs/sleep_frailty_ml_neuronal_activity'


def integer_coordinate(value):
    d=Decimal(value)
    if not d.is_finite() or d!=d.to_integral_value():raise ValueError(value)
    return int(d)


def main():
    extra={}
    with gzip.open(OUT/'original_neuronal_activity_identity_inventory.csv.gz','rt') as handle:
        for r in csv.DictReader(handle):
            if r['history_class']=='new_physical_sequence':extra[r['name']]=r
    counts=collections.Counter();categories=collections.Counter();classes=collections.Counter()
    sources=collections.Counter();lengths=collections.Counter();builds=collections.Counter();unique=set()
    canonical={'chr'+str(i) for i in range(1,23)}|{'chrX','chrY'}
    with gzip.open(OUT/'ngn2_element_activity_metadata_202602_bbmap_only_mpralib095_normalized_bcalm_dna_ge_1.tsv.gz','rt') as handle:
        for r in csv.DictReader(handle,delimiter='\t'):
            if r['name'] not in extra:continue
            categories[r['category']]+=1;classes[r['class']]+=1;sources[r['source']]+=1
            lengths[len(r['sequence'])]+=1;builds[r['ref']]+=1
            unique.add(extra[r['name']]['physical_sequence_sha256_both_orientations'])
            try:
                a,b=integer_coordinate(r['start']),integer_coordinate(r['end'])
                valid=r['chr'] in canonical and a>=0 and b>a
            except (ValueError,InvalidOperation):valid=False
            counts['has_declared_canonical_genomic_interval' if valid else 'no_declared_canonical_genomic_interval']+=1
            if valid and r['ref']=='GRCh38':counts['declared_GRCh38_canonical_interval']+=1
    receipt={'new_literal_source_construct_rows':len(extra),'new_literal_distinct_full_DNA_sequences':len(unique),
       'category':dict(categories),'class':dict(classes),'source':dict(sources),'lengths':dict(lengths),'declared_reference':dict(builds),
       'coordinate_availability':dict(counts),'coordinate_parser':'Exact Decimal integral parsing accepts 2181843.0; never truncates a nonintegral coordinate.',
       'biologically_independent_new_regions_established':False,'new_validation_data':False,'current_training_changed':False,
       'remaining':'Genomic component and neighborhood exclusions across every development source and the whole reserved study; duplicated construct treatment; raw-count/historical-software limitations remain. No new target is inserted into running fits.'}
    (OUT/'additional_activity_construct_categories.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt),flush=True)


if __name__=='__main__':main()
