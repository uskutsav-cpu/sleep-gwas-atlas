#!/usr/bin/env python3
"""Create additive source-only PRELAUNCH metadata artifacts; no genotype hashing/decoding."""
import csv
import datetime
import hashlib
import json
import stat
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
P = ROOT / 'sleep_unified_research_v4'
S = P / 'source_provenance'
GENO = Path('/Volumes/Extreme SSD/brain6-work/brain6-alternative-local-validation-v1/reference')
NOW = datetime.datetime.now(datetime.timezone.utc).isoformat()


def regular(p):
    p = Path(p)
    assert not p.is_symlink() and all(not q.is_symlink() for q in p.parents), str(p)
    assert stat.S_ISREG(p.stat().st_mode), str(p)
    return p


def small_sha(p):
    p = regular(p)
    assert p.stat().st_size <= 1048576 and p != GENO / '1000G_EUR.bed'
    return hashlib.sha256(p.read_bytes()).hexdigest()


def small_json(p):
    small_sha(p)
    return json.loads(Path(p).read_text())


def ident(p, header_bytes=0, small_content=False):
    p = regular(p)
    st = p.stat()
    result = dict(path=str(p), bytes=st.st_size, device=st.st_dev, inode=st.st_ino, mtime_ns=st.st_mtime_ns)
    if header_bytes:
        with p.open('rb') as f:
            header = f.read(header_bytes)
        result.update(header_read_bytes=len(header), header_sha256=hashlib.sha256(header).hexdigest())
        if p.suffix == '.bed':
            assert header_bytes == 3
            result['header_hex'] = header.hex()
    result['fresh_full_content_sha256'] = small_sha(p) if small_content else None
    return result


def save(p, obj):
    payload = (json.dumps(obj, indent=2) + '\n').encode()
    expected = hashlib.sha256(payload).hexdigest()
    with p.open('xb') as f:
        f.write(payload)
    assert small_sha(p) == expected


inputs = [
    P / 'reviews/genomicsem_realistic_ld_calibration_feasibility_v4.md',
    P / 'manifests/canonical_realistic_ld_source_admission_plan_v4.json',
    P / 'statistical_validation/canonical_realistic_ld_feasibility_probe_receipt_v4.json',
    P / 'statistical_validation/COMMON_BOUNDARY_METHOD_VALIDATION_PROTOCOL_v4.md',
    P / 'manifests/canonical_200_interval_preparation_v4.json',
    P / 'statistical_validation/canonical_200_intervals_v4.tsv',
    P / 'scripts/41_prepare_canonical_blocks.py',
    P / 'statistical_validation/canonical_signed_genotype_generator_fixture_v4.py',
    P / 'statistical_validation/canonical_signed_genotype_generator_fixture_receipt_v4.json',
    P / 'statistical_validation/canonical_ld_feasibility_sources_v4/zenodo_8292725_metadata.json',
    P / 'statistical_validation/canonical_ld_feasibility_sources_v4/integrated_call_samples_v3.20130502.ALL.panel',
    ROOT / 'config/interpretation_ldsc_seg_reference.json',
    P / 'manifests/global_SSD_resource_reservation_v4_4.json',
    P / 'manifests/global_SSD_resource_reservation_v4_5.json',
    P / 'logs/closed_sensitivity_storage_reservation_release_v4_5.json',
    GENO / 'prepared_v1/reference_prep_provenance.json',
    GENO / 'per_chr_v1/split_reference_provenance.json',
    GENO.parent / 'SUPERGNOVA_source/README.md',
    Path('/Users/swethasunilkumar/Documents/Codex/2026-10-08/goa/work/sleep-gwas-atlas/brain6/scripts/prepare_supergnova_reference_v1.py'),
    Path('/Users/swethasunilkumar/Documents/Codex/2026-10-08/goa/work/sleep-gwas-atlas/brain6/results/brain6_alternative_local_validation_v1/reference_prep_stdout.json'),
]
before = {str(p): small_sha(p) for p in inputs}
assert before[str(P / 'manifests/global_SSD_resource_reservation_v4_5.json')] == '7ce69095a6f67ae681d547476d1557afac2de0f97d083b12bc1b689d0788417f'
assert before[str(P / 'logs/closed_sensitivity_storage_reservation_release_v4_5.json')] == 'ec0362b6edb45d7cb7d4cde701d9e6a1c27c6eaa684d3743171e57a9e4b28566'
assert before[str(P / 'statistical_validation/COMMON_BOUNDARY_METHOD_VALIDATION_PROTOCOL_v4.md')] == '4cc651a2b51b63cf9ea1e7535da48960c7d395967879e8bfe2fb781855542445'
prior_plan = small_json(P / 'manifests/canonical_realistic_ld_source_admission_plan_v4.json')
prior_probe = small_json(P / 'statistical_validation/canonical_realistic_ld_feasibility_probe_receipt_v4.json')
local_prep = small_json(GENO / 'prepared_v1/reference_prep_provenance.json')
ledger = small_json(P / 'manifests/global_SSD_resource_reservation_v4_5.json')
assert ledger['reserved_total_bytes'] == 319399338749 < ledger['ceiling_bytes'] == 322122547200
assert ledger['component_bytes']['proposed_signed_reference_source_verification_bytes'] == 3221225472
assert ledger['component_bytes']['sensitivity_bytes'] == 268435456
assert ledger['closed_sensitivity_reexecution_admitted'] is False
assert prior_plan['restricted_S_qualification_admitted'] is False and prior_plan['full_C_calibration_admitted'] is False
bed = ident(GENO / '1000G_EUR.bed', 3)
bim = ident(GENO / '1000G_EUR.bim', 4096)
fam = ident(GENO / '1000G_EUR.fam', small_content=True)
afreq = ident(GENO / '1000G_EUR.afreq', 256)
assert bed['header_hex'] == '6c1b01' and bed['bytes'] == 231387159
assert fam['fresh_full_content_sha256'] == local_prep['source_sha256']['1000G_EUR.fam']
assert all(bed[k] == prior_probe['raw_bed_metadata_only'][k] for k in ['bytes','device','inode','mtime_ns'])
fam_rows = [line.split() for line in (GENO / '1000G_EUR.fam').read_text().splitlines()]
assert all(len(r) == 6 for r in fam_rows) and len(fam_rows) == 503
assert len({r[1] for r in fam_rows}) == 503
panel_path = P / 'statistical_validation/canonical_ld_feasibility_sources_v4/integrated_call_samples_v3.20130502.ALL.panel'
with panel_path.open(newline='') as f:
    panel = list(csv.DictReader(f, delimiter='\t'))
assert {r[1] for r in fam_rows} == {r['sample'] for r in panel if r['super_pop'] == 'EUR'}
zenodo_path = P / 'statistical_validation/canonical_ld_feasibility_sources_v4/zenodo_8292725_metadata.json'
zenodo = small_json(zenodo_path)
zf = next(f for f in zenodo['files'] if f['key'] == '1000G_Phase3_plinkfiles.tgz')
assert zf['size'] == 288277344 and zf['checksum'] == 'md5:a7773ab485827b533cb300c76356d76b'
assert zenodo['metadata']['license']['id'] == 'cc-by-4.0'
namespace = ledger['new_source_verification_namespace']
qc = {
    'scope': 'REFERENCE_ONLY_NO_GWAS_Z_P_RG_OR_NATIVE_FINAL_MASK_SELECTION',
    'source_precondition': 'AUTHENTICATED_UPSTREAM_ARCHIVE_OR_EXACT_REPRODUCIBLE_ORIGINAL_CONTENT_CHAIN_REQUIRED_FOR_SOURCE_ADMISSION',
    'raw_source_only_not_prepared_MAF05_or_cm_interpolated_copies': True,
    'reference_universe': '22 exact consumed LD-score ID/CHR/BP rows and exact coordinate/allele map; no score/weight/M substitution',
    'master_rank': 'Store rank0 and rank1 for every historical master reference row before MHC/static/source-QC exclusion; chromosome numeric1..22 and original row order. C03 future cycling uses original rank0, never subset renumbering; simulation admission remains false.',
    'coordinate_map': {'path':prior_plan['coordinate_map_path'], 'sha256':prior_plan['coordinate_map_sha256']},
    'reference_sha256': prior_plan['reference_sha256'],
    'primary_interval': {'path':prior_plan['primary_interval_path'], 'sha256':prior_plan['primary_interval_sha256']},
    'MHC': 'Exclude chr6 BP25000000 through34000000 inclusive, exactly the existing canonical builder; retain master ranks and exclusion records.',
    'static_eligibility': 'Unique ID with exact CHR/BP and exact ordered or swapped A1/A2; both alleles distinct one-base ACGT; exclude A/T,C/G; no complement-strand rescue; missing allele-map/ID, coordinate/allele conflict and duplicate source IDs recorded/rejected.',
    'dosage_A1': {'00':2, '01':'MISSING', '10':1, '11':0},
    'samples': 'Same503 in source FAM order; all22 FAM files byte-identical and authoritative Phase3 EUR sample-set match; do not drop/reorder samples to pass QC.',
    'padding': 'For503 samples every final variant byte high2 padding bits must be zero; wrong mode/size/truncation/trailing genotype bytes is fatal.',
    'variant_missingness': 'Retain only missing_count<=25 out503; >25/503 excluded before universe seal. All-missing is separately classified and excluded.',
    'zero_variance': 'Exclude if observed_nonmissing_count<2 or exact integer statistic nobs*sum_dose2-sum_dose**2<=0; no floor/jitter/clipping.',
    'MAF': 'Observed A1 frequency and MAF logged; no additional MAF threshold and no AFREQ-based selection. No HWE/outcome/power filter.',
    'missing_imputation': 'Impute each retained SNP to observed mean; center all503; scale by sqrt(sum(centered**2)/503). Require finite values, zero column mean and unit Gram diagonal.',
    'allele_sign': 'Map dosage-counted BIM A1 to target coordinate-map A1: +1 ordered match, -1 swapped match. Apply sign to centered standardized column; record both ordered allele pairs and sign. No assumption that A1 is minor/reference/effect allele.',
    'source_order': 'Record chromosome/member/BIM row0 and BED byte offset3+row0*ceil(503/4), full master rank0/rank1 and final J index; preserve ties and deterministic CHR/BP/ID order.',
    'frozen_universe_output': 'Ordered exact J and sample-order digests; every static and genotype-QC exclusion with reason; never repair a failed source by an outcome-informed subset.',
}
controls = {
    'tiny_format_controls': 'Hand-encoded all4 genotype states,503-sample final-byte padding, truncated/trailing/sample-major inputs, wrong FAM order, A1/A2 swap, mean missingness, exact zero/all-missing and near-zero nonzero variance; discrete counts/alleles/offsets exact.',
    'decoder_independence': 'Before actual-source PASS, compare bit-shift decoderA to separately maintained pinned decoderB (e.g. independently reviewed PLINK allele-preserving export). Do not describe two wrappers around one decoder as independent. Executable/runtime/version/code hash and A1 export orientation still require freeze; PATH census found no plink/plink2.',
    'source_windows': 'After reference-only J freeze: for each chromosome select first64, middle64 and last64 eligible J rows, no overlaps; fewer192 is source-control INCOMPLETE. Also use fixed32/32/32 early/middle/late scatter to retain distant within-chromosome correlations; no favorable window substitution.',
    'direct_operator': 'For each selected source window, independently form float64 G and R=G.T@G/503; compare deterministic-vector direct Rv to two-pass G.T@(Gv/503), including distant rows, chunk sizes1/7/64/2048 and reversed chunk iteration with unchanged SNP/sample order.',
    'sign_control': 'Flip predetermined every-third selected target allele; require Gnew=G*D, Rnew=D*R*D and same sample-centered norm; squared-r agreement alone is insufficient.',
    'global_operator': 'Streaming full-chromosome deterministic probes, energy v.T R v=||Gv||**2/503 and bilinear symmetry; no dense global SNP matrix and no hard jackknife/storage-window boundary truncation.',
    'tolerance': 'Fixed protocol point/deletion tolerance abs_error<=1e-10+1e-8*abs(reference); discrete decode/count/sign checks exact; unit-column mean/norm absolute1e-10; no numerical repair to pass.',
    'PSD': 'Gram factor is PSD by construction; within-chromosome centered rank<=502. Dense-window numeric checks do not establish population LD/long-range adequacy.',
    'controls_executed_here': 0,
}
plan = {
    'schema':'signed_ld_source_authentication_QC_prelaunch_plan_v4_1', 'prepared_utc':NOW,
    'admitted':False, 'execution_authorized':False, 'body_worker_implemented_or_reviewed':False,
    'source_only_metadata_checker_code':str(Path(__file__).resolve()),
    'scope':'SOURCE_AUTHENTICATION_AND_REFERENCE_ONLY_QC_PRELAUNCH_DESIGN_NO_SIMULATION_OR_ESTIMATOR',
    'current_local_source': {'BED':bed,'BIM':bim,'FAM':fam,'AFREQ':afreq,'historical_expected_BED_sha256':local_prep['source_sha256']['1000G_EUR.bed'],'raw_variant_count_from_prior_manifest':1836406,'upstream_genotype_content_authenticated':False,'fresh_BED_or_BIM_full_hash':False},
    'unresolved_original_source_chain':'No exact raw-genotype transfer URL/release/checksum/archive-member/extraction receipt or reproducible official source-to-compact-subset derivation was found in the bounded inspected metadata. Local prep/split hashes and FAM match do not supply it.',
    'separate_authenticated_candidate':{'record_id':8292725,'record_doi':'10.5281/zenodo.8292725','creator':'Steven Gazal','metadata_path':str(zenodo_path),'metadata_sha256':small_sha(zenodo_path),'file_id':zf['id'],'file_name':zf['key'],'content_url':zf['links']['self'],'bytes':zf['size'],'upstream_md5':zf['checksum'],'observed_record_license':'CC BY4.0','record_is_exact_version_not_latest_concept':True,'payload_acquired':False,'archive_member_identity_not_yet_known':True,'no_EAS_GRCh38_or_mirror_substitution':True,'fresh_local_SHA256_required_before_and_after_use':True},
    'archive_members': {'regular_file_count':66,'chromosomes':list(range(1,23)),'basenames':'1000G.EUR.QC.{1..22}.{bed,bim,fam}','allowed_uniform_prefixes':['','1000G_EUR_Phase3_plink/'],'prefix_qualification':'Archive packaging unexamined: require exactly one uniform frozen prefix family; unexpected names/entries stop rather than silently changing scope.','directories':'Only root and matching single approved prefix; reject other directories.','reject':['duplicate names','absolute/traversal paths','symlinks','hardlinks','devices','FIFOs','sparse members','unbounded extended metadata','unexpected regular files','unverified trailing gzip payload'],'extraction':'Manual owned-target regular file writes only; no extractall, owner/mode/xattr preservation or overwrite. Bind every member size/SHA and archive before/after SHA/MD5.'},
    'strict_source_checks':['All22 FAM files byte-identical and503 unique IIDs/503 authoritative Phase3 EUR samples; preserve exact order and FID/IID digest.','Each BED6c1b01 and exact3+nvariant*126bytes; no truncation/trailing/padding violation.','Each BIM6fields, expected numeric chromosome, positive 1-based BP, nondecreasing position; source IDs globally unique or reject before choosing J.','Build cannot be proved by directory name: document primary Phase3/GRCh37 evidence and exact selected ID/CHR/BP comparison to verified GRCh37 map; any conflict is retained exclusion; unresolved full-source build remains qualified, not asserted.','No assumption new archive equals old compact BED or historical eur_w_ld_chr; prove exact content derivation separately if claiming equivalence.'],
    'resource_contract':{'ledger_path':str(P/'manifests/global_SSD_resource_reservation_v4_5.json'),'ledger_sha256':small_sha(P/'manifests/global_SSD_resource_reservation_v4_5.json'),'global_ceiling_bytes':322122547200,'reserved_total_bytes':319399338749,'unallocated_margin_bytes':2723208451,'source_namespace':namespace,'source_aggregate_cap_bytes':3221225472,'archive_and_all_failed_partial_cap_bytes':536870912,'extracted_all_partial_and_failed_cap_bytes':2147483648,'QC_controls_metadata_and_all_failed_output_cap_bytes':536870912,'canonical16control8GiB_reused':False,'closed_sensitivity_cap_bytes':268435456,'closed_sensitivity_reexecution_admitted':False,'worker_count':1,'aggregate_RSS_stop_bytes':2147483648,'BLAS_threads':1,'decoded_chunk_max_SNPs':2048,'internal_floor_bytes':3221225472,'SSD_floor_bytes':5368709120,'shared_heavy_mutex':prior_plan['future_resource_contract_review_required']['required_heavy_lock'],'whole_stage_monotonic_deadline_seconds':7200,'current_whole_campaign_byte_meter_required_at_every_stage':True,'meter_includes_regular_files_AppleDouble_failed_partials_and_pending_outputs':True,'runtime_copy_installation':False,'output_paths_SSD_only':True,'lock_policy':'Existing exact shared heavy mutex nonblocking exclusive for owned worker lifetime; do not create a separate lock or wait indefinitely.','guard_policy':'Reviewed controller checks RSS/floors/deadline/component/subcaps/current300GiB meter and CLOSEDsens<=256MiB every1second during owned child/hash/extraction/QC and before persistence; stop only owned process/group, defer catchable signals during scoped cleanup, preserve failure/partial receipt. Current scanner/meter/runtime/controller code still requires independent freeze/review.'},
    'minimal_executable_future_stage_design':[
        {'stage':'G0','operation':'Freeze actual controller/runtime/decoder hashes and exact reviewed source plan/admission receipt. Recheck all dependencies, regular leaves/all parents, component/global current census, CLOSEDsens bindings/size and existing floors; acquire shared heavy mutex.','failure':'NO_SOURCE_EXECUTION_ADMISSION'},
        {'stage':'G1','operation':'Optional separately admitted local-continuity check:64KiB streaming hash raw BED/BIM/FAM before/after plus stat/header. Compare expected local hashes.','qualification':'Even exact fresh96da… proves continuity only, not upstream content authentication.'},
        {'stage':'G2','operation':'If newly admitted: one TLS-verified HTTPS-only bounded owned curl fetch of exact record8292725 EUR archive, no automatic retry/resume/reuse; bind fresh official metadata/HEAD size, full exact288277344bytes/upstreamMD5/newSHA before publication to owned cache.','failure':'Retain failed partial under512MiB, no automatic replacement archive.'},
        {'stage':'G3','operation':'Guarded bounded tar/gzip header/member pass and manual66-file extraction to new owned subdirectory; exact allowlist; prove member hashes, gzip completion, FAM/BED/BIM strict checks and archive unchanged.','failure':'SOURCE_ARCHIVE_OR_MEMBER_IDENTITY_REJECTED'},
        {'stage':'G4','operation':'Freeze static reference/coordinate intersections and whole master ranks, decode source-only QC in chunks<=2048, save every exclusion and exact J/sample/row/sign digest under512MiBQCcap.','failure':'SOURCE_QC_REJECTED_OR_INCOMPLETE_NO_UNIVERSE_SUBSTITUTION'},
        {'stage':'G5','operation':'Execute frozen independent decoder/direct-window/sign/global-operator controls; persist fixed intended immutable source/QC seal with full input/output identity checks before/after each terminal write.','failure':'SOURCE_CONTROL_REJECTED_OR_INCOMPLETE'},
    ],
    'reference_only_QC':qc, 'independent_controls':controls,
    'finite_reference_truth_contract':{'R':'G.T@G/503 within chromosome, zero between chromosomes; fixed503 samples/mean imputation/centering/scaling/allele signs; no block boundary factorization or eigenvalue clipping.','rank':'Each chromosome rank<=502; sample-induced long-range correlations remain present within chromosome.','meaning':'Conditional-on-authenticated-finite-reference Gaussian model, not250000 participant genotypes or population-LD validation.','future_integrated_truth':'If S later admitted, integrated random-effect Sigma_g ratios, signal covarianceR^2 tensorSigma_g/m_total and noiseR tensorSigma_e; no realized-effect normalization/target switch. Variant-N scaling does not establish real cohort-overlap laws.','historical_scores_weights_M':'All exact historical files retained; no reconstructed replacements. Finite R^2 need not match their L2/window/bias-correction interpretation.','not_established':['actual pair retained-set coverage','population LD','physical variant-N overlap','liability transforms','winner selection','full C calibration']},
    'preserved_full_protocol':prior_plan['preserved_protocol'],
    'quota_or_acceptance_margin_reduction':False,'S_simulations_run':0,'C_calibration_credit':0,'estimator_calls':0,'biological41Cov_release_admitted':False,
    'source_stage_launch_prerequisites_remaining':['Independent actual acquisition/extraction/hash/QC/controller code review and fixed intent seals','Frozen existing runtime and independently maintained decoderB identity','Root explicit source-only admission after exact resource/current-meter preflight','Upstream/archive content authentication and selected-build compatibility before source admission'],
}
save(P/'manifests/signed_ld_503eur_source_only_prelaunch_plan_v4_1.json',plan)
evidence_rows=[]
def row(id,proof,meaning,gap,url=''):
    evidence_rows.append({'evidence_id':id,'metadata_proof':proof,'established':meaning,'not_established_or_gate':gap,'primary_url':url,'observed_utc':NOW})
row('LOCAL_BED',json.dumps(bed,sort_keys=True),'Exact current path/stat/three-byte mode/dimension metadata agrees prior probe; historical expectedSHA96da… exists.','No fresh full hash, decode or upstream content authentication.')
row('LOCAL_BIM',json.dumps(bim,sort_keys=True),'Current small header and stat; prior verified fullSHA0375…/1836406rows consumed as metadata.','No fresh full BIM hash/read; coordinates/alleles not independently re-established for whole source here.')
row('LOCAL_FAM',json.dumps(fam,sort_keys=True),'Fresh smallFAM SHAccc5…;503uniqueIDs equal authoritative entirePhase3EURsample set.','FAM ancestry does not authenticate BED genotypes or upstream conversion.')
row('LOCAL_PREPARATION',str(GENO/'prepared_v1/reference_prep_provenance.json')+' | SHA256='+small_sha(GENO/'prepared_v1/reference_prep_provenance.json'),'Local raw/prepared input/output hashes and historical MAF/nonpalindromic/map filtering recorded.','Original source URL/checksum/archive/member/extraction/derivation absent; source genotype hash is local, not publisher authentication.')
row('LOCAL_CHROMOSOME_SPLIT',str(GENO/'per_chr_v1/split_reference_provenance.json')+' | SHA256='+small_sha(GENO/'per_chr_v1/split_reference_provenance.json'),'Downstream exact-copy claim/member digests recorded.','Descendant records cannot independently authenticate unbound ancestor.')
row('SUPERGNOVA_SOURCE_README',str(GENO.parent/'SUPERGNOVA_source/README.md')+' | SHA256='+small_sha(GENO.parent/'SUPERGNOVA_source/README.md'),'Author README generically points to1000G bfiles FTP and says rare variants filtered.','No exact raw231387159byte file/member checksum or release binding; pinned README web fetch failed.','https://github.com/qlu-lab/SUPERGNOVA/blob/319e84e114a4f954005a4592756c56cfee083667/README.md')
row('PHASE3_PANEL',str(panel_path)+' | SHA256='+small_sha(panel_path),'Authoritative sample labels support503EUR FAM identities only.','No genotype content equivalence or conversion proof.','https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel')
row('ZENODO_EXACT_CANDIDATE',str(zenodo_path)+' | SHA256='+small_sha(zenodo_path),'Official record8292725 EUR file288277344B MD5a7773…;creator/DOI/licenseCCBY4.0.','Archive/member/genotypes not acquired or verified; no proof of identity with old compact BED.','https://zenodo.org/records/8292725')
row('DISTINCT_MIRROR_CLAIM',str(ROOT/'config/interpretation_ldsc_seg_reference.json')+' | SHA256='+small_sha(ROOT/'config/interpretation_ldsc_seg_reference.json'),'Pinned songlab mirror config lists66files/1593837597B.','Distinct source metadata claim; not actual Zenodo content/member identity proof.')
row('GLOBAL_LEDGER_V4_5',str(P/'manifests/global_SSD_resource_reservation_v4_5.json')+' | SHA256='+small_sha(P/'manifests/global_SSD_resource_reservation_v4_5.json'),'Separate3GiB source cap fits unchanged300GiB;CLOSEDsens256MiB release record bound.','No source or worker admission; current whole-campaign meter required at each execution boundary.')
row('PLINK_FORMAT','Dated official web facts; no HTML body SHA invented.','Published two-bit/BIM/FAM order defines planned A1 dosage/padding controls.','DecoderB/current binary not frozen; no decoder/source control executed.','https://www.cog-genomics.org/plink/1.9/formats#bed')
with (S/'signed_ld_503eur_provenance_evidence_v4_1.tsv').open('x',newline='') as f:
    writer=csv.DictWriter(f,fieldnames=list(evidence_rows[0]),delimiter='\t');writer.writeheader();writer.writerows(evidence_rows)
after={str(p):small_sha(p) for p in inputs};assert before==after
assert ident(GENO/'1000G_EUR.bed',3)==bed and ident(GENO/'1000G_EUR.bim',4096)==bim
assert ident(GENO/'1000G_EUR.fam',small_content=True)==fam
receipt={'schema':'signed_ld_503eur_small_metadata_header_receipt_v4_1','recorded_utc':NOW,'script_sha256':small_sha(__file__),'input_metadata_sha256_before':before,'input_metadata_sha256_after':after,'raw_source_header_and_stat':{'BED':bed,'BIM':bim,'FAM':fam,'AFREQ':afreq},'BED_bytes_read_per_probe':3,'BIM_bytes_read_per_probe':4096,'AFREQ_bytes_read_per_probe':256,'fresh_full_BED_or_BIM_SHA_or_genotype_decode':False,'official_panel_FAM_EUR_set_equal':True,'local_upstream_genotype_content_authenticated':False,'new_archive_downloaded_or_read':False,'new_source_QC_controls_or_simulations_executed':0,'estimator_calls':0,'protected_project_writes':0,'full_calibration_credit':0,'new_global3GiB_allocation_is_not_execution_admission':True,'web_primary_sources_checked':['https://www.cog-genomics.org/plink/1.9/formats#bed','https://zenodo.org/records/8292725','https://github.com/bulik/ldsc/wiki/LD-Score-Estimation-Tutorial'],'web_access_failures':['PinnedSUPERGNOVA README GitHub/Raw cache miss; local metadata README retained.'],'full_protocol_read_and_8000_quota_preserved':True,'outputs_exclusive':True}
save(S/'signed_ld_503eur_metadata_header_receipt_v4_1.json',receipt)
print(json.dumps({'plan':str(P/'manifests/signed_ld_503eur_source_only_prelaunch_plan_v4_1.json'),'plan_sha256':small_sha(P/'manifests/signed_ld_503eur_source_only_prelaunch_plan_v4_1.json'),'receipt_sha256':small_sha(S/'signed_ld_503eur_metadata_header_receipt_v4_1.json'),'upstream_authenticated':False,'launch_admitted':False}))
