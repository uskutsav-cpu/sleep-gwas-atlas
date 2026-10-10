#!/usr/bin/env python3
"""Independent small forensic checks; never fits GWAS or runs calibration.

Reads pinned code, compact archived numerical tables, and bounded reference
headers. GWAS source integrity is inherited explicitly from the v4 execution
plan, not re-certified here. Writes only new specialist v4 filenames.
"""
from pathlib import Path
from datetime import datetime, timezone
import csv
import gzip
import hashlib
import json
import math

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'sleep_unified_research_v4/statistical_validation'
SRC = OUT / 'genomicsem_specialist_sources_v4'
ARCHIVE = Path('/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/FAILED-PARTIAL-2026-09-01-sleep_gwas_atlas')
read_hashes = {}


def read(path):
    path = Path(path)
    data = path.read_bytes()
    read_hashes[str(path)] = hashlib.sha256(data).hexdigest()
    return data


def table(path):
    return list(csv.DictReader(read(path).decode().splitlines(), delimiter='\t'))


def cov_manual(x, y):
    assert len(x) == len(y) and len(x) > 1
    mx, my = sum(x) / len(x), sum(y) / len(y)
    return sum((a - mx) * (b - my) for a, b in zip(x, y)) / (len(x) - 1)


def ratio(a, c, b):
    assert a > 0 and b > 0
    return c / math.sqrt(a * b)


def partial_scalar(a, b, c, u, v, w):
    assert w > 0
    ar, br, cr = a - u*u/w, b - v*v/w, c - u*v/w
    assert ar > 0 and br > 0
    return cr / math.sqrt(ar * br)


def r_index_blocks(coords, count):
    # Exact mathematical equivalent of v0.1.1's one-based floor(seq(...)).
    n = len(coords)
    starts = [math.floor(1 + (n-1)*i/count) for i in range(count+1)]
    ends = [s-1 for s in starts[1:count]] + [n]
    return [coords[a-1:b] for a, b in zip(starts[:count], ends)]


def main():
    plan = json.loads(read(ROOT/'sleep_unified_research_v4/manifests/ssd_native_execution_plan_v4_3.json'))
    pairs = table(ARCHIVE/'results/tables/ldsc_covariance_pairs.tsv')
    diagonal = {r['trait_1']: r for r in pairs if r['trait_1'] == r['trait_2']}
    fixed_errors = []
    for r in pairs:
        den = math.sqrt(float(diagonal[r['trait_1']]['genetic_covariance']) * float(diagonal[r['trait_2']]['genetic_covariance']))
        expected = float(r['genetic_covariance_se'])/den
        fixed_errors.append(abs(float(r['genetic_correlation_se'])-expected)/expected)
    model_fit = table(ARCHIVE/'results/tables/genomic_sem_model_fit.tsv')
    model_provenance = json.loads(read(ARCHIVE/'results/tables/factor_gwas.provenance.json'))
    diagnostics = table(ARCHIVE/'results/tables/ldsc_covariance_diagnostics.tsv')
    metadata = table(ARCHIVE/'results/tables/ldsc_covariance_metadata.tsv')
    map_path=ARCHIVE/'work/track_b_completion/local_dependency_copies/ref/hm3_grch37_variant_map.tsv.gz'
    map_provenance=json.loads(read(map_path.with_name(map_path.name+'.provenance.json')))
    map_hash=hashlib.sha256()
    with map_path.open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''):map_hash.update(chunk)
    read_hashes[str(map_path)]=map_hash.hexdigest()
    assert map_hash.hexdigest()==map_provenance['map_sha256']
    with gzip.open(map_path,'rt') as f:map_header=next(f).strip().split('\t')
    matched_refs=all(plan['dependencies_sha256'].get(str(map_path.parent/'eur_w_ld_chr'/r['path']))==r['sha256']
                     for r in map_provenance['source']['files'] if r['path'].endswith('.l2.ldscore.gz'))
    assert matched_refs
    for p in [ROOT/'scripts/25_genomicsem_covariance.R', ROOT/'scripts/28_chromosome_split_covariance.R',
              ROOT/'scripts/29_genomicsem_model.R', ROOT/'scripts/30_finalize_genomicsem.py',
              ROOT/'sleep_unified_research_v1/scripts/native_ldsc_capture.py',
              ROOT/'sleep_unified_research_v1/SHARED_ESTIMATOR_COVARIANCE.md',
              ROOT/'sleep_unified_research_v1/reviews/statistical_validity_v1.json',
              ROOT/'sleep_unified_research_v4/reviews/statistical_geneticist_critical_review_v4.md',
              ROOT/'config/liftover_plans.tsv',ROOT/'environment/tool_versions.tsv',
              ROOT.parent/'ldsc-code/ldscore/jackknife.py', ROOT.parent/'ldsc-code/ldscore/irwls.py',
              ROOT.parent/'ldsc-code/ldscore/regressions.py', ROOT.parent/'ldsc-code/ldscore/sumstats.py']:
        read(p)
    for p in SRC.iterdir():
        if p.is_file(): read(p)
    refs = sorted((Path(p),h) for p,h in plan['dependencies_sha256'].items() if p.endswith('.l2.ldscore.gz'))
    headers = []
    for p,h in refs:
        with gzip.open(p,'rt') as f:
            header = next(f).strip().split('\t')
            first = next(f).strip().split('\t')
        assert set(('CHR','SNP','BP')).issubset(header)
        headers.append({'path':str(p), 'build':'GRCh37/hg19_pinned_reference_contract',
                        'expected_sha256':h, 'hash_verification_scope':'inherited_fresh_v4_3_plan_not_rehashed_here',
                        'header':'|'.join(header), 'first_chr':first[header.index('CHR')],
                        'first_snp':first[header.index('SNP')], 'first_bp':first[header.index('BP')],
                        'all_rows_coordinate_uniqueness':'NOT_YET_SCANNED_IN_THIS_REVIEW',
                        'common_boundary_manifest':'NOT_YET_CREATED'})
    with (OUT/'genomicsem_block_coordinate_availability_v4.tsv').open('w') as f:
        w = csv.DictWriter(f,fieldnames=list(headers[0]),delimiter='\t');w.writeheader();w.writerows(headers)

    # Synthetic examples are independent of every GWAS outcome.
    coords1=list(range(1,801)); coords2=[x for x in coords1 if not 101<=x<=200]
    blocks1=r_index_blocks(coords1,4);blocks2=r_index_blocks(coords2,4)
    assert any(a!=b for a,b in zip(blocks1,blocks2))
    # Shared coordinate intervals preserve deletion meaning despite missing SNPs.
    canonical=[(1,201),(201,401),(401,601),(601,801)]
    joined1=[[x for x in coords1 if lo<=x<hi] for lo,hi in canonical]
    joined2=[[x for x in coords2 if lo<=x<hi] for lo,hi in canonical]
    assert all(set(b).issubset(a) for a,b in zip(map(set,joined1),map(set,joined2)))

    a,c,b=.3,.09,.2
    rv=ratio(a,c,b)
    analytic=[-rv/(2*a),1/math.sqrt(a*b),-rv/(2*b)]
    vals=[a,c,b]; fd=[]
    for i in range(3):
        plus=vals.copy();minus=vals.copy();eps=1e-6
        plus[i]+=eps;minus[i]-=eps
        fd.append((ratio(*plus)-ratio(*minus))/(2*eps))
    gradient_error=max(abs(x-y) for x,y in zip(analytic,fd))
    assert gradient_error<1e-8
    assert ratio(.3,.3,.3)==1

    # Two independent scalar derivations: covariance residualisation vs
    # the standardized one-confounder partial-correlation identity.
    a,b,c,u,v,w=.3,.2,.09,.06,.04,.25
    p1=partial_scalar(a,b,c,u,v,w)
    rab=ratio(a,c,b);rac=ratio(a,u,w);rbc=ratio(b,v,w)
    p2=(rab-rac*rbc)/math.sqrt((1-rac*rac)*(1-rbc*rbc))
    assert abs(p1-p2)<1e-14
    # Unknown/negative residual variance must fail, never be clipped or smoothed.
    rejected=False
    try: partial_scalar(.01,.2,.09,.06,.04,.25)
    except AssertionError: rejected=True
    assert rejected

    delete_a=[.22,.20,.23,.19];delete_b=[.20,.21,.20,.19]
    n=len(delete_a);full_a=.21;full_b=.20
    pseudo_a=[n*full_a-(n-1)*x for x in delete_a]
    pseudo_b=[n*full_b-(n-1)*x for x in delete_b]
    ca=cov_manual(pseudo_a,pseudo_a)/n; cb=cov_manual(pseudo_b,pseudo_b)/n
    cab=cov_manual(pseudo_a,pseudo_b)/n
    vc=ca+cb-2*cab
    direct=cov_manual([x-y for x,y in zip(pseudo_a,pseudo_b)],[x-y for x,y in zip(pseudo_a,pseudo_b)])/n
    deletion_formula=(n-1)**2/n*cov_manual(delete_a,delete_b)
    assert abs(vc-direct)<1e-15 and abs(cab-deletion_formula)<1e-15
    assert abs(cab)<=math.sqrt(ca*cb)+1e-15

    report={'utc':datetime.now(timezone.utc).isoformat(),
        'scope':'Independent forensic arithmetic/code review; no native fit, empirical calibration, heterogeneity outcome, conditional GWAS result, S/V smoothing or manuscript',
        'status':'COMMON_BOUNDARY_METHOD_IMPLEMENTABLE_BUT_NOT_CALIBRATED_OR_ADMITTED_FOR_INFERENCE',
        'native_campaign_inputs':{'count':len(plan['inputs_verified']), 'hash_match_count':sum(x['match'] for x in plan['inputs_verified']),
            'verification_scope':'inherited_v4_3_plan; only bounded reference headers reread, not 158 complete GWAS files'},
        'reference_headers_with_exact_chr_snp_bp':len(headers),
        'existing_coordinate_map':{'path':str(map_path),'bytes':map_path.stat().st_size,
            'sha256':map_hash.hexdigest(),'matches_archived_provenance':True,
            'all22_map_source_reference_hashes_match_campaign':matched_refs,
            'header':map_header,'genome_build':map_provenance['genome_build'],
            'mapped_hm3_ids_provenance_count':map_provenance['mapped_rows'],
            'unmapped_hm3_ids_provenance_count':map_provenance['hm3_rows_without_eur_ldscore_coordinate'],
            'row_count_scope':'provenance count; compressed bytes freshly verified, only header reread',
            'actual_final_pair_snp_join_verified':False},
        'archived_export_arithmetic':{'rows':len(pairs),'fixed_denominator_se_max_relative_error':max(fixed_errors),
            'diagonal_rg_identically_one':all(float(r['genetic_correlation'])==1 for r in diagonal.values()),
            'diagonal_nonzero_exported_rg_se':sum(float(r['genetic_correlation_se'])>0 for r in diagonal.values())},
        'archived_diagnostics_retrieved_not_reeigendecomposed':diagnostics,
        'archived_metadata':metadata,
        'archived_heldout_models':{'count':len(model_fit), 'validated_count':sum(r['validation_status']=='VALIDATED' for r in model_fit),
            'negative_residual_model_count':sum(int(r['negative_observed_residuals'])>0 for r in model_fit),
            'terminal_status':model_provenance['terminal_status']},
        'small_independent_checks':{'pair_specific_index_blocks_misaligned':True,
            'example_index_boundary_chr1_bp_end_A':[x[-1] for x in blocks1],
            'example_index_boundary_chr1_bp_end_B':[x[-1] for x in blocks2],
            'canonical_coordinate_intervals_respect_missingness':True,
            'ratio_gradient_finite_difference_max_abs_error':gradient_error,
            'partial_correlation_two_formula_max_abs_error':abs(p1-p2),
            'nonpositive_residual_variance_rejected':rejected,
            'contrast_variance_identity_max_abs_error':abs(vc-direct),
            'pseudovalue_delete_covariance_identity_max_abs_error':abs(cab-deletion_formula)},
        'primary_pins':{'partialLDSC':'v0.1.1 @ 18b743baf609ac3898b00ee92f37783dc3b61347',
            'GenomicSEM':'0.0.5 @ 6b65ca5db39fdade08b0d811477be1cdd57b5039',
            'CBIIT_LDSC':'6c673952cee74bd5c57aef1555a03b1c015399a0'},
        'empirical_calibration_achieved':False, 'joint_covariance_estimated':False,
        'common_boundary_manifest_created':False,'new_primary_hypothesis_admitted':False,
        'read_input_sha256':read_hashes}
    report['script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (OUT/'genomicsem_specialist_forensic_checks_v4.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:report[k] for k in ['status','reference_headers_with_exact_chr_snp_bp','archived_export_arithmetic','archived_heldout_models','small_independent_checks']},indent=2))


if __name__=='__main__': main()
