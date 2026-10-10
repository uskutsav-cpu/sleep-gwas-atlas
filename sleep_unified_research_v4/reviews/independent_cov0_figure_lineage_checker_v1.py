"""Additive stdlib figure lineage/SVG text-and-geometry audit; no plotting rerun.

Reads sealed review metadata, one aggregate table, the producer source, manifest,
root visual QA and six figure exports. Does not repeat41-row statistical
equations/native audits or read interpreter/package/GWAS/reference bodies.
"""
import ast
import csv
import hashlib
import json
from pathlib import Path
import re
import resource
import struct
import textwrap
import time
import xml.etree.ElementTree as ET

P=Path(__file__).resolve().parents[1]
R=P/'reviews'
TABLE=P/'tables/shared_sleep_cov0_diagnostic41_v4_1.tsv'
PRODUCER=P/'scripts/107_export_cov0_diagnostic_figures_v4_1.py'
MANIFEST=P/'manifests/shared_sleep_cov0_diagnostic_figure_exports_v4_1.json'
QA=P/'logs/shared_sleep_cov0_diagnostic_visual_QA_v4_1.json'
OLD_SEAL=R/'independent_shared_sleep_cov0_diagnostic41_review_seal_v1.json'
OUT=R/'independent_cov0_figure_lineage_receipt_v1.json'
started=time.monotonic()
fixed={
    str(TABLE):'0698b9493637176c03b975bc36e9fbd6e2c5baf6dfdc545486c004df693149bb',
    str(PRODUCER):'b5f224b1121d09cde24c4aa157a2216e562eff2d732d051c00089ddbc1a2628d',
    str(MANIFEST):'33e9450d04d62483bc286b440e3d32e226729b300ce288b759f5f953e14396c9',
    str(QA):'92787d6a4a1b29a780c6ca4e6dea3b55c5a3f1c3612039bc27ace82f6646fff0',
    str(OLD_SEAL):'fac58ba84ae14620672bb1c236ab6b3caaedb1b6e4094213b77fa8f3ac11317b',
}


def sha(path):
    path=Path(path)
    assert path.is_file() and not path.is_symlink()
    h=hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda:stream.read(65536),b''):h.update(chunk)
    return h.hexdigest()


def bound():
    assert time.monotonic()-started<60
    assert resource.getrusage(resource.RUSAGE_SELF).ru_maxrss<(64 << 20)


assert all(sha(q)==v for q,v in fixed.items())
old=json.loads(OLD_SEAL.read_text())
assert old['file_sha256'][str(TABLE)]==fixed[str(TABLE)]
assert old['verdict']=='QUALIFIED_AGGREGATE_ARITHMETIC_AND_JOIN_PASS'
# Only hash metadata receipt/report files from the old arithmetic seal.
# No equation checker or old test suite is invoked.
assert len(old['file_sha256'])==13
assert all(sha(q)==v for q,v in old['file_sha256'].items())
with TABLE.open(newline='') as stream:rows=list(csv.DictReader(stream,delimiter='\t'))
assert len(rows)==len({r['pair_id'] for r in rows})==41
manifest=json.loads(MANIFEST.read_text())
qa=json.loads(QA.read_text())
assert manifest['script_sha256']==fixed[str(PRODUCER)] and manifest['table_sha256']==fixed[str(TABLE)]
assert manifest['figure_page_count']==2 and manifest['export_count']==6
assert manifest['sampling_covariance_established'] is False and manifest['confirmed_heterogeneity'] is False
assert manifest['new_BH_or_other_correction'] is False and manifest['new_estimator_fits']==0
assert manifest['manuscript_legends_authored'] is False
assert manifest['pycache_writes_disabled'] is True
assert manifest['visual_QA_pending'] is True  # Retained export-time state.
assert qa['manifest_sha256']==fixed[str(MANIFEST)] and qa['table_sha256']==fixed[str(TABLE)]
assert qa['root_viewed_page_count']==2 and qa['visual_QA_pending'] is False
assert qa['sampling_covariance_or_heterogeneity_validated'] is False
assert qa['arithmetic_independent_review_is_separate'] is True
assert qa['new_primary_claim_or_manuscript_authored'] is False
for key in ['readable_labels_and_numerical_intervals','no_clipping_or_overlap_seen',
            'same_axis_scale_both_pages','unverified_Cov0_and_reused_sleep_qualification_on_each_page']:
    assert qa[key] is True
ids=[r['pair_id'] for r in rows]
assert qa['exact41_pair_ids']==ids
assert [f['page'] for f in manifest['figures']]==[1,2]
assert [f['supported_value_count'] for f in manifest['figures']]==[21,20]
assert [x for f in manifest['figures'] for x in f['pair_ids']]==ids
exports={q:v for f in manifest['figures'] for q,v in f['output_sha256'].items()}
assert len(exports)==6 and exports==qa['all6_export_sha256']
assert all(sha(q)==v for q,v in exports.items())
source=PRODUCER.read_text()
tree=ast.parse(source)
labels=next(ast.literal_eval(n.value) for n in tree.body if isinstance(n,ast.Assign)
            and any(isinstance(t,ast.Name) and t.id=='LABELS' for t in n.targets))
main=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='main')
strings=[n.value for n in ast.walk(main) if isinstance(n,ast.Constant) and isinstance(n.value,str)]
assert not any(s.startswith('difference_P') or s.startswith('difference_Z') or s.startswith('current_cov0_nominal') or s.startswith('original_nominal7') for s in strings)
assert 'difference_external_minus_discovery' in strings
assert 'difference_CI95_low_cov0_unvalidated' in strings and 'difference_CI95_high_cov0_unvalidated' in strings
assert "forest.set_xlim(*limits)" in source
assert source.count("color='#263f58'")==2
assert "range(0, len(rows), 21)" in source and "selected = rows[start:start+21]" in source
low=min(float(r['difference_CI95_low_cov0_unvalidated']) for r in rows)
high=max(float(r['difference_CI95_high_cov0_unvalidated']) for r in rows)
span=max(high-low,0.5)
limits=(min(low-.06*span,-.03),max(high+.06*span,.03))
assert limits[0]<0<limits[1]
NS={'s':'http://www.w3.org/2000/svg'}
captions=[
    'Discovery–validation effect differences',
    'Unverified Cov=0 diagnostic; the same sleep GWAS is reused',
    'All 41 estimated locked pairs retained. Shared sampling covariance is unestimated; these intervals cannot establish heterogeneity.',
    'Table: shared_sleep_cov0_diagnostic41_v4_1.tsv | No new test family, primary claim or fully independent two-trait replication.',
]
pages=[]
max_geometry_error=0.0
for page,f in enumerate(manifest['figures'],1):
    bound()
    selected=rows[(page-1)*21:page*21]
    assert f['pair_ids']==[r['pair_id'] for r in selected]
    assert f['scientific_scope']=='UNVERIFIED_ZERO_COVARIANCE_DIAGNOSTIC_ONLY'
    assert f['visual_QA_pending'] is True and f['all_text_inside_canvas'] is True
    assert f['SVG_title_description_role'] is True and f['PDF_tagged_accessibility_claimed'] is False
    assert f['minimum_text_points']==9.5
    assert len(f['output_sha256'])==3
    bysuffix={Path(q).suffix:Path(q) for q in f['output_sha256']}
    assert set(bysuffix)=={'.pdf','.svg','.png'}
    root=ET.parse(bysuffix['.svg']).getroot()
    assert root.attrib['role']=='img'
    assert root.attrib['aria-labelledby']==f'shared_sleep_cov0_diagnostic_page{page}-title shared_sleep_cov0_diagnostic_page{page}-desc'
    title=root.find('s:title',NS).text
    desc=root.find('s:desc',NS).text
    assert title=='Qualified zero-covariance effect-difference diagnostic'
    assert f'Page {page} shows {len(selected)} of all41' in desc
    for snippet in ['unverified Cov=0 standard error','The sleep GWAS is reused.',
                    'No covariance-adjusted inference or confirmed heterogeneity is supported.']:
        assert snippet in desc
    alltext=[''.join(x.itertext()) for x in root.findall('.//s:text',NS)]
    assert all(s in alltext for s in captions)
    assert f'Page {page}/2' in alltext
    assert 'External rg − discovery rg; interval assumes Cov=0' in alltext
    axes={x.attrib.get('id'):x for x in root.findall('.//s:g',NS)}
    left,numbers=axes['axes_1'],axes['axes_3']
    lefttext=[''.join(x.itertext()) for x in left.findall('.//s:text',NS)]
    expectedleft=['Sleep trait | outcome']
    for r in selected:
        expectedleft+=textwrap.fill(labels[r['sleep_trait']]+' | '+r['phenotype'],width=45,
            break_long_words=False,break_on_hyphens=False).splitlines()
    assert lefttext==expectedleft
    expectednumbers=['Difference [95% CI0]']+[f"{float(r['difference_external_minus_discovery']):.3f} [{float(r['difference_CI95_low_cov0_unvalidated']):.3f}, {float(r['difference_CI95_high_cov0_unvalidated']):.3f}]" for r in selected]
    numbertext=[''.join(x.itertext()) for x in numbers.findall('.//s:text',NS)]
    assert numbertext==expectednumbers
    forest=axes['axes_2']
    patch=next(g for g in forest.findall('./s:g',NS) if g.attrib.get('id')=='patch_2')
    coords=[float(x) for x in re.findall(r'-?\d+(?:\.\d+)?',patch.find('s:path',NS).attrib['d'])]
    xleft,xright=coords[0],coords[2]
    assert (xleft,xright)==(418.32,776.16)
    ci=[];dots=[]
    for g in forest.findall('./s:g',NS):
        if not g.attrib.get('id','').startswith('line2d_'):continue
        for path in g.findall('./s:path',NS):
            if 'stroke: #263f58' in path.attrib.get('style',''):
                values=[float(x) for x in re.findall(r'-?\d+(?:\.\d+)?',path.attrib['d'])]
                assert len(values)==4 and values[1]==values[3]
                ci.append(values)
        for use in g.findall('.//s:use',NS):
            if 'stroke: #263f58' in use.attrib.get('style',''):
                assert 'fill: #263f58' in use.attrib.get('style','')
                dots.append((float(use.attrib['x']),float(use.attrib['y'])))
    assert len(ci)==len(dots)==len(selected)
    geom=[]
    def project(value):return xleft+(value-limits[0])/(limits[1]-limits[0])*(xright-xleft)
    for rowline,(r,bar,dot) in enumerate(zip(selected,ci,dots),(page-1)*21+2):
        errors=[abs(bar[0]-project(float(r['difference_CI95_low_cov0_unvalidated']))),
                abs(bar[2]-project(float(r['difference_CI95_high_cov0_unvalidated']))),
                abs(dot[0]-project(float(r['difference_external_minus_discovery']))),abs(dot[1]-bar[1])]
        assert max(errors)<1e-6  # SVG coordinate serialization only, in points.
        max_geometry_error=max(max_geometry_error,max(errors))
        geom.append(dict(pair_id=r['pair_id'],table_row=rowline,interval_x=[bar[0],bar[2]],
            dot_x=dot[0],max_geometry_error_points=max(errors)))
    with bysuffix['.png'].open('rb') as stream:header=stream.read(24)
    assert header[:8]==b'\x89PNG\r\n\x1a\n'
    dimensions=list(struct.unpack('>II',header[16:24]))
    assert dimensions==f['PNG_dimensions_pixels']
    assert all(abs(d-300)<.02 for d in f['PNG_dpi'])
    with bysuffix['.pdf'].open('rb') as stream:assert stream.read(5)==b'%PDF-'
    pages.append(dict(page=page,members=len(selected),exact_ordered_pair_ids=f['pair_ids'],
        SVG_labels_numbers_qualifiers_geometry_pass=True,geometry=geom,
        common_x_limits=list(limits),PNG_dimensions=dimensions,
        independent_raster_or_PDF_visual_inspection=False,root_visual_QA_separate=True))
assert all(sha(q)==v for q,v in fixed.items()) and all(sha(q)==v for q,v in exports.items())
bound()
result=dict(schema='independent_cov0_figure_lineage_metadata_svg_review_v1',all_checks_pass=True,
    fixed_current_inputs_sha256=fixed,export_sha256=exports,pages=pages,total_members=41,
    page_sizes=[21,20],same_x_limits=list(limits),max_geometry_error_points=max_geometry_error,
    geometry_tolerance_points=1e-6,geometry_tolerance_is_export_serialization_not_scientific=True,
    exact_plotted_fields=['difference_external_minus_discovery','difference_CI95_low_cov0_unvalidated','difference_CI95_high_cov0_unvalidated'],
    plotted_nominal_P_classification_or_historical7_significance_colors=False,
    new_significance_primary_or_independent_replication_claim=False,
    qualifiers_present_in_actual_each_SVG_text_and_description=True,
    arithmetic_seal_unchanged=True,old13_metadata_artifacts_current=True,
    arithmetic_equations_328_native190_or_old_suite_reruns=0,
    GWAS_reference_nativevector_runtime_body_reads=False,network_locks_workers_fits_censuses=0,
    root_visual_QA_identity_current=True,independent_raster_review=False,
    current_recorded_plotting_environment=manifest['environment'],
    environment_versions_binary_or_cache_use_independently_attested=False,
    prior_exports_or_manifest_pending_flags_modified=False,
    new_inference_sampling_covariance_or_confirmed_heterogeneity_validated=False,
    elapsed_seconds=time.monotonic()-started,max_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    limitations='Source/metadata/SVG coordinate/text lineage plus current export hashes and distinct root visual QA. PDF/PNG hashed and signature/dimensions checked; no independent raster/PDF visual review. Same plotting source and final export identities support lineage, not native/input/covariance or package-binary certification.')
with OUT.open('x') as stream:json.dump(result,stream,indent=2,allow_nan=False);stream.write('\n')
print(json.dumps(dict(all_pass=True,pages=2,members=41,exports=6,x_limits=list(limits),
    max_geometry_error_points=max_geometry_error,elapsed_seconds=result['elapsed_seconds'],max_RSS_bytes=result['max_RSS_bytes'])))
