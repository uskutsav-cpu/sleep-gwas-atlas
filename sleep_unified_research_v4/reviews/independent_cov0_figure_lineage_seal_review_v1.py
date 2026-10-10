"""Seal additive exact-export lineage review; no plot or statistical rerun."""
import datetime
import hashlib
import json
from pathlib import Path

P=Path(__file__).resolve().parents[1]
R=P/'reviews'
PREFIX='independent_cov0_figure_lineage'
RECEIPT=R/(PREFIX+'_receipt_v1.json')
CHECKER=R/(PREFIX+'_checker_v1.py')
REPORT=R/(PREFIX+'_review_v1.json')
MD=R/(PREFIX+'_review_v1.md')
SEAL=R/(PREFIX+'_review_seal_v1.json')


def sha(path):
    path=Path(path)
    assert path.is_file() and not path.is_symlink()
    h=hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda:stream.read(65536),b''):h.update(chunk)
    return h.hexdigest()


def save(path,value):
    with path.open('x') as stream:json.dump(value,stream,indent=2,allow_nan=False);stream.write('\n')


r=json.loads(RECEIPT.read_text())
assert r['all_checks_pass'] is True and r['total_members']==41 and r['page_sizes']==[21,20]
files=dict(r['fixed_current_inputs_sha256'])
files.update(r['export_sha256'])
assert all(sha(q)==v for q,v in files.items())
oldpath=R/'independent_shared_sleep_cov0_diagnostic41_review_seal_v1.json'
old=json.loads(oldpath.read_text())
assert sha(oldpath)=='fac58ba84ae14620672bb1c236ab6b3caaedb1b6e4094213b77fa8f3ac11317b'
assert all(sha(q)==v for q,v in old['file_sha256'].items())
report=dict(schema='independent_cov0_figure_lineage_additive_review_v1',
    reviewed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
    verdict='QUALIFIED_FIGURE_LINEAGE_PASS',material_lineage_or_qualification_errors=[],
    exact_current_inputs_sha256=r['fixed_current_inputs_sha256'],six_current_exports_sha256=r['export_sha256'],
    additive_only_preserved_arithmetic_seal_sha256='fac58ba84ae14620672bb1c236ab6b3caaedb1b6e4094213b77fa8f3ac11317b',
    exact41_complete_ordered_membership=True,page_sizes=[21,20],common_x_limits=r['same_x_limits'],
    source_and_actual_SVG_text_geometry_lineage_pass=True,
    exact_plotted_fields=r['exact_plotted_fields'],
    all41_actual_SVG_markers_and_CI_segments_checked=True,
    max_SVG_coordinate_error_points=r['max_geometry_error_points'],serialization_tolerance_points=1e-6,
    numeric_labels_round_only_to3decimals=True,
    significance_colors_or_P_labels_primary_promotion=False,
    original_classification_and_historical7_remain_in_linked_table_not_new_plot_selection=True,
    each_page_unverified_Cov0_reused_sleep_covariance_unestimated_no_heterogeneity_primary_family_independence_qualifiers=True,
    root_visual_QA_bound_and_distinct=True,independent_raster_or_PDF_visual_QA=False,
    PDF_and_PNG_hashes_and_signatures_checked=True,actual_PNG_dimensions_checked=True,
    SVG_title_description_role_verified=True,PDF_tagged_accessibility_claim=False,
    recorded_plotting_environment=r['current_recorded_plotting_environment'],
    runtime_binary_packages_cache_behavior_independently_verified=False,
    old_exports_or_pending_flags_overwritten=False,
    completed_statistical_equations_native190_old_suites_repeated=False,
    biological_reference_nativevector_runtime_body_reads=False,
    network_native_worker_fit_mutex_or_resource_census_calls=0,
    sampling_covariance_estimated_confirmed_heterogeneity_source_admission_or_new_scientific_claim=False,
    resources=dict(elapsed_seconds=r['elapsed_seconds'],max_RSS_bytes=r['max_RSS_bytes']),
    limitations=r['limitations'])
save(REPORT,report)
producer=P/'scripts/107_export_cov0_diagnostic_figures_v4_1.py'
manifest=P/'manifests/shared_sleep_cov0_diagnostic_figure_exports_v4_1.json'
qa=P/'logs/shared_sleep_cov0_diagnostic_visual_QA_v4_1.json'
table=P/'tables/shared_sleep_cov0_diagnostic41_v4_1.tsv'
text=f"""# Additive all41 Cov0 figure lineage review

**QUALIFIED_FIGURE_LINEAGE_PASS.** No source/metadata/annotation/geometry lineage error was found in the exact [107 producer]({producer}), [manifest]({manifest}) and six exports. The completed arithmetic [seal]({oldpath}) (`fac58ba84ae14620672bb1c236ab6b3caaedb1b6e4094213b77fa8f3ac11317b`) and its13 artifacts remain unchanged. This additive check does not repeat328 decimal equations, any native190 audit or old suite.

The manifest SHA is `33e9450d04d62483bc286b440e3d32e226729b300ce288b759f5f953e14396c9`, producer `b5f224b1121d09cde24c4aa157a2216e562eff2d732d051c00089ddbc1a2628d`, and distinct root [visual QA]({qa}) `92787d6a4a1b29a780c6ca4e6dea3b55c5a3f1c3612039bc27ace82f6646fff0`. All six PDF/SVG/PNG export hashes currently match both the manifest and QA. The [independent checker]({CHECKER}) and [receipt]({RECEIPT}) bind these identities to the exact41-row [table]({table}), SHA `0698b9493637176c03b975bc36e9fbd6e2c5baf6dfdc545486c004df693149bb`.

The pages contain21 and20 pairs, respectively, with every pair ID in exact table order and no omissions, duplication, filtering or ranking by P. The plotting source selects only `difference_external_minus_discovery` and the two `difference_CI95_*_cov0_unvalidated` fields. Each actual SVG was parsed independently using the standard library. All41 ordered labels and3-decimal numeric annotations agree with the table; all41 dot coordinates and41 interval segments agree with the common table-based linear axis transform within1e-6 points, a vector-coordinate serialization allowance. Maximum coordinate error is4.99×10^-7 points. Both forest panels use the same x limits, [-0.40733258369902675,0.326278068463656], and the same fixed dot/interval color. This is a check of figure coordinates, not re-estimation of the statistical values or a scientific tolerance.

There are no significance colors, P labels, stars or primary-result promotion. The figures show all41 differences/CI0; the original23/18 classification and seven-row annotation remain in the linked table rather than being used for a new plot selection. Each SVG's actual text states the same sleep GWAS is reused, Cov=0 is unverified, shared covariance is unestimated, the intervals cannot establish heterogeneity, and no new test family, primary claim or fully independent two-trait replication follows. Titles/descriptions, role and accessible labeling agree with that scope. PDF tagged accessibility is explicitly not claimed.

PDF and PNG files were hashed and their signatures checked; actual PNG dimensions agree with the manifest. I did not independently raster-review the PDF/PNG pages. Root's bound QA separately records viewing both pages, readable labels/numerical intervals, no clipping or overlap and a common scale. Source/SVG lineage and this recorded visual QA remain distinct evidence. Export-time `visual_QA_pending=True` fields are preserved; the additive root QA records its later `False` state without rewriting the earlier manifest.

The current manifest records Python3.13.5, Matplotlib3.11.1, NumPy2.5.3 and Pillow12.3.0. This review preserves those recorded versions and root's read-only -B/chat-cache qualification without reading runtime/package bodies or independently attesting the binary versions/cache behavior. Older exports are untouched. No GWAS/reference/native-vector bodies, fits, workers, network, heavy mutex or resource/runtime census were used; the bounded checker completed in0.023seconds with24,854,528 bytes peak RSS.

These pages are supported as explicitly qualified diagnostic displays of the existing marginal-input table. The lineage/visual evidence does not estimate sampling covariance, validate heterogeneity, clear source/clinical semantics, establish independent replication, or admit a new biological, causal, novelty or primary scientific claim. The statistical limitations in the preserved arithmetic review continue to apply.
"""
with MD.open('x') as stream:stream.write(text)
files.update(old['file_sha256'])
for q in [CHECKER,RECEIPT,REPORT,MD,Path(__file__)]:files[str(q)]=sha(q)
assert all(sha(q)==v for q,v in files.items())
save(SEAL,dict(schema='independent_cov0_figure_lineage_additive_seal_v1',
    verdict=report['verdict'],file_sha256=files,regular_file_count=len(files),
    preserves_prior_arithmetic_seal=True,independent_root_visual_QA_distinct=True,
    no_completed_equations_native190_or_old_suite_repeat=True,source_covariance_or_new_biological_inference_admitted=False))
assert all(sha(q)==v for q,v in files.items())
print(json.dumps(dict(verdict=report['verdict'],md=str(MD),md_sha256=sha(MD),
    report=str(REPORT),report_sha256=sha(REPORT),seal=str(SEAL),seal_sha256=sha(SEAL),file_count=len(files))))
