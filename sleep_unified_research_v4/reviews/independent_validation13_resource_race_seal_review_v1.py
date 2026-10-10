"""Seal targeted race diagnosis and result-free continuation requirements."""
import datetime
import hashlib
import json
from pathlib import Path

P=Path(__file__).resolve().parents[1]
R=P/'reviews'
PREFIX='independent_validation13_resource_race'
CHECKER=R/(PREFIX+'_checker_v1.py')
RECEIPT=R/(PREFIX+'_receipt_v1.json')
MD=R/(PREFIX+'_review_v1.md')
REPORT=R/(PREFIX+'_review_v1.json')
SPEC=R/(PREFIX+'_continuation_requirements_v1.json')
SEAL=R/(PREFIX+'_review_seal_v1.json')
ROOT_STOP=P/'logs/validation_raw13_v3_actual_stop_root_receipt_v1.json'
ROOT_STOP_SHA='e1142ad40f71b73c0af1e0cd2526902628b614d74df04561b2ed46a1261fd606'


def sha(path):
    path=Path(path)
    assert path.is_file() and not path.is_symlink()
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path,value):
    with path.open('x') as stream:json.dump(value,stream,indent=2,allow_nan=False);stream.write('\n')


r=json.loads(RECEIPT.read_text())
assert r['measured_operational_race_reproduced_with_exact_meter_AST'] is True and r['control_count']==5
assert sha(ROOT_STOP)==ROOT_STOP_SHA
root=json.loads(ROOT_STOP.read_text())
assert root['source_bodies_completed']==11 and root['whole_family_terminal_success'] is False
assert root['pending_retained'] is True and root['source13_started'] is False
assert root['failed_partial_bytes']==57671680 and root['failed_partial_sha256']==r['source12_prefix_sha256_recorded_not_freshly_rehashed']
for entry in r['completed11_source_receipt_bindings_from_failed_family_only']:
    assert root['file_sha256'][entry['path']]==entry['sha256']
# Intentionally do not traverse root's30 referenced files/body-prefix binding.
files=dict(r['exact_failed_source_family_and_resource_code_metadata_sha256'])
files[str(ROOT_STOP)]=ROOT_STOP_SHA
assert all(sha(q)==v for q,v in files.items())
spec=dict(schema='targeted_validation13_additive_continuation_requirements_v1',
    status='RESULT_FREE_REQUIREMENTS_ONLY_NOT_IMPLEMENTATION_OR_ADMISSION',
    old_epoch_preserve='v3 ownership/PENDING/failedfamily/two source12 failurecopies/logs/partial/header bindings; do not overwrite/remove/convert tosuccess.',
    old_plan_sha256=r['recorded_v3_plan_sha256'],old_admission_sha256=r['recorded_admission_sha256'],
    old_failure_metadata_sha256=files,
    future_namespace='New distinctly versioned validation_raw_replay_v4 namespace and frozen executor/plan/admission/receipts; not prepared or executed bythisreview.',
    unchanged_scientific_source_membership='13 historical sources in exact original order, URLs/generations/fullbodySHA/MD5/bytes/sourcequalifications;217/41 statistical scopes unchanged.',
    inherited_successful11_immutable_receipt_binding=r['completed11_source_receipt_bindings_from_failed_family_only'],
    old11_not_reread_or_current_body_certified_in_this_review=True,
    old11_adoption_conditions=[
        'Explicit new root admission permits individual authenticated source adoption despite old familyFAILED/PENDING; does not declare oldwholefamilycomplete.',
        'Exact old successful receiptSHA beforeparse/afterconsumption, both copies, originalmember/plan/source identity, full originalSHA_MD5_size, HTTP generation/ETag/length, return0/stopnull/ownedreapedteardown/failureaddendumveto.',
        'Fresh current sourcebody fullSHA_MD5_size by separately guarded root/newexecutor verification; no copy/move/hardlink or newdownload of11oldbodies.',
        'Freeze per-source originplan/receipt/bodypath identity in newall13master; no pretending adopted11receiptshave newepochplanSHA.',
    ],
    new_downloads_only=dict(indices=[12,13],commands=2,source12_original_member=r['source12_expected_original_member'],
        source13='Unchanged exact member13 from frozen original13plan/design; not newlyinspected here.',
        mode='FULL_HTTP200_NEW_BODY_ONLY',Range_resume_admitted=False,automatic_retry=False),
    meter_correction_contract=[
        'New version only; within exact owned campaign root, require root exists and isregular directory with trustedownership/no symlink traversal.',
        'Single no-follow metadata stat perentry; count observed regular st_size, skip symlinks, traverse dirs without following replaced symlinks.',
        'Catch only descendantENOENT from concurrently disappearing file/directory during observational census; log disappearedpath/count and observedbyte total. Root disappearance remainsfatal.',
        'PermissionError/EIO/othererrors, malformed sizes, actual capoverrun and interrupted/resourcefailed scans remainfatal; no broadexcept=>zero/stale-success fallback.',
        'This is non-atomic resourceobservation, not contentidentity certification or proof of instantaneous global occupancy. Keep frozen reservation/componentcaps and all floors/poll/finalguards.',
        'Do not catch or forgive missingchangedbound sources/receipts/outputs/plan/protocol/body prefix/currentfullbody hashes; strictimmutable evidence callbacks retainfailures.',
    ],
    proposed_resource_delta=dict(exact_failedsource12_prefix_bytes=57671680,
        proposed_fixed_new_operational_metadata_cap_bytes=16<<20,total_proposed_allocation_bytes=57671680+(16<<20),
        proposal_requires_author_check_and_freeze=True,
        allocate_from='Prior4_5unallocatedmargin in separatelyversioned newreservationledger only; unchanged300GiB ceiling, oldactiveledger and activeworkers untouched.',
        original13_declared_body_bytes=10058648185,original13_plus_this_prefix_minimum_bytes=10116319865,
        original13_plus_prefix_and_proposed_metadata_bytes=10133097081,
        no_currentcampaign_or_margin_census_made_here=True,no_duplicate11bodyreservation=True),
    unchanged_limits=r['unchanged_required_caps'],
    deadlines='Explicitly freeze the continuation deadline and old elapsedtime/start; recommend carrying forward original96h familywall deadline rather than silently resetting it. Each new fullbody retains2h worker ceiling.',
    ownership='Same exclusive singletransfer familyflock; deferredSIGINT/TERM/HUP; ownedcurlgroup terminated/reaped/quarantined beforeunlock; allowedonecurl alongsideexistingoneheavyworker only under unchangedfrozenresourcebudget.',
    final_success_oracle='All13 exactcurrent authenticatedmembers in originalorder:11 adopted+2new, immutable per-source originclosure, two identical primary receipts with precomputedintendedSHA, Terminal2privatePENDING/exactseal, currentsource/header/hash/resource/deferredsignal gates before/afterpersistence; no full13success beforeall13.',
    prelaunch_test_requirements=[
        'Healthy stablecensus; disappearance beforestat/afterobservedstat/withinchilddir doesnotabortsafeownedcensus; symlink exclusion; capoverrun and permission/EIO/rootloss stillreject.',
        'Strictbound-evidence ENOENT and mutation stillreject; no softened source/receipt/output identitygate.',
        'Adopt11wrongstatus/command/member/originSHA/bodySHA/HTTP/cleanup/failureaddendum casesreject; no absentoldPENDING pretense or oldfamilypromotion.',
        'Exact two fullHTTP200 commands, no Range/append/retry; end-to-end fulloriginalSHA_MD5_bytes for12/13; all13 cardinality/order/oracle verified.',
        'Precomputed two-primary/sealintent, postpersistence evidence/resource/signal faults, failedsupplementalwrites, ownedcleanup ordering and nofalsewholefamilycommit.',
    ],
    author_and_reviewer='Different authoragent prepares replacement; this reviewer independently reviews exactnewepochplan/code beforeseparaterootadmission.',
    scientific_source_QC_power_replication_admission=False)
save(SPEC,spec)
report=dict(schema='independent_validation13_targeted_resource_race_review_v1',
    reviewed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
    verdict='RESOURCE_CENSUS_ENOENT_RACE_CONFIRMED_CONTINUATION_NOT_ADMITTED',
    diagnosed_gap='Exact90 globalcensus doubleobserves is_file andstat, allowing concurrentephemeraldescendantdeletion tothrowFileNotFoundError.',
    independent_exact_original_meter_controls=5,independent_originalmeter_race_witness=True,
    current_source_family_failed_evidence_sha256=files,
    root_stop_receipt_bound=True,root_partial_and30metadataverification_is_distinct_from_ourtargetedreads=True,
    eleven_complete_sources_recorded_not_freshly_body_reverified=True,
    preserved_source12_prefix_bytes=57671680,preserved_prefix_SHA_from_receipts=r['source12_prefix_sha256_recorded_not_freshly_rehashed'],
    curl_returncode=-15,emptygroup_ownedteardown_verified_inreceipt=True,
    recorded_floors_healthy=True,actual_globalcap_atstop_notmeasured=True,
    source13_notstarted_inrootreceipt=True,
    current_v3_completedfull13_or_pipeline_scientific_success=False,
    safe_additive_continuation_spec_path=str(SPEC),safe_additive_continuation_spec_sha256=sha(SPEC),
    continuation_full12and13_only=True,Range_resume_admitted=False,new_source_substitution=False,
    originalfailedfamily_PENDING_overwrite_or_automaticretry_admitted=False,
    newexecutor_orcommonmeter_orfutureplan_implemented=False,separatefutureprelaunchreview_androotadmission_required=True,
    no11body_or_successful11receipt_reread_no_realcampaigncensus_orworker_transfer_oldsuite=True,
    scientific_QC_or_source_clearance_or_replication_claim=False,
    resource_observations=dict(elapsed_seconds=r['elapsed_seconds'],max_RSS_bytes=r['max_RSS_bytes']))
save(REPORT,report)
code=P/'scripts/90_acquire_validation_raw_sources_v3.py'
text=f"""# Independent targeted validation13 resource-race diagnosis

**RESOURCE_CENSUS_ENOENT_RACE_CONFIRMED; continuation is not yet admitted.** The actual stoppedv3 metadata and exact [90 acquisition resource-meter function]({code}:88) support an observational census race. The meter evaluates `p.is_file()` and then a separate `p.stat().st_size` while another permitted worker removes an ephemeral descendant. A disappeared path between these observations propagates FileNotFoundError. The transfer then correctly fail-stops and preserves evidence. This diagnosis does not show a GWAS-source identity/QC error or recertify completed science.

The two actual family receipt copies are identical, SHA `995a00139dfbb40f772dfdd262ad93f20fc9fbef65acddf9b2cb50d46e3a6a3d`; they record11 completed source receipts and `FAMILY_STOP_PRESERVED_REQUIRES_EXPLICIT_AUDIT`. The exact recorded plan is `8e93d548163bda284b6292b9780962060bbc3ff4afde44e5a916c80f3b69b708`. The two source12 K11_CHOLELITH receipt copies are identical, SHA `bde5ece0dfb1aa799502a7f57ca57b6856e7c0130fd7ccbf3964784ec73ae662`. They record FileNotFoundError on `core_pipeline/large35_bounded_replay_v3/napping/ephemeral_bounded_spool/columns/typed_08_00450000.pkl`, a6.153second stop, curl20956 SIGTERM/returncode-15, owned teardown verified and empty remaining group. The inspected monitor cleanup waits for the owned process; this is receipt/code evidence, not a fresh OS process census.

Source12's failed prefix is57,671,680 bytes, recorded SHA `3537a1b14c7e3560c43d38898dfadd6ba280738498ef5de4c21195c962a1a415`. Its exact expected full source is809,815,233 bytes, SHA `50af445a7edcd53129da6e4f21a461af62ae9dd53ad983f00525bb2f29a298e1`, MD5 `29db68b6a610c94cb5b8d48b7b01ccb7`, pinned GCS generation1777990004571514. I did not reread the prefix, any11 full bodies or their successful receipts. Root's separate [stop receipt]({ROOT_STOP}), SHA `{ROOT_STOP_SHA}`, reports fresh30 metadata/prefix bindings and that source13 never started; this current root record is bound without traversing its referenced body files.

The actual source/family post-stop snapshots record approximately17.03GB internal free space and741GB SSD free space, above the unchanged3GiB/5GiB floors; peak observed curl RSS was7,602,176 bytes and deadlines were not reached. The census itself failed, so this review does not assert a measured global300GiB occupancy at the failure instant. Failed/PENDING/ownership records and originals remain intact; no current worker or active ledger was modified.

The [independent checker]({CHECKER}) extracts only the exact90 meter AST and runs five bounded metadata controls: a stable tiny regular file is counted; a private one-byte file removed after is_file produces the same FileNotFoundError; PermissionError remains fatal; an over-cap count rejects; an excluded symlink does not contribute. The checker traverses a private wrapper rather than the live SSD campaign and does not invoke acquisition, old suites or corrected production code. The [control receipt]({RECEIPT}) and all tiny witnesses are sealed. The measured run used0.445seconds and26,296,320 bytes RSS.

The correction must be additive and scoped to observational resource counting in the trusted owned campaign. Use one no-follow metadata stat for each observed entry, retain its regular-file size, and handle only descendantENOENT from disappearing entries/directories, recording disappearance diagnostics. Root loss, permissions/EIO/other errors and actual cap overruns remain fatal. Do not follow symlink replacements or broadly substitute zero/stale totals on failure. A concurrent census is not an atomic snapshot; the frozen reservation/component caps, floors and polling/final guards remain required. Missing or changed bound inputs, receipts, plans, outputs, source hashes and prefix identity must still fail their strict evidence consumers. No production patch or corrected-meter admission is made by this report.

The result-free [continuation specification]({SPEC}) freezes the safe direction agreed by root: a distinctly versioned epoch, no Range or resume, only new fullHTTP200 acquisitions for source12 and source13, and explicit authenticated reuse of the existing11 originals without copy/move/new download. New admission must qualify adopting individual exact successful source proofs from the failed old family without claiming that old family completed. Freeze each original receipt SHA, plan/member/endpoint/body path and successful HTTP/cleanup oracle, then separately verify current complete SHA/MD5/bytes under guarded resources. Preserve the old11 receipt origin identities rather than rewriting them to the new plan. The new master must retain all13 original source identities/order once, with11 inherited and2 new complete authenticated inputs.

Budget the exact57,671,680-byte failed prefix plus a fixed modest new operational metadata allowance from the4_5 unallocated margin in a new reviewed ledger, under the unchanged300GiB ceiling. The specification proposes16MiB metadata (74,448,896 bytes total additional allocation), subject to author verification/freezing; it does not assume new unlimited capacity or independently attest the remaining margin. The original13 declared body total remains10,058,648,185 bytes; with this retained prefix it is at least10,116,319,865 bytes, before the new metadata allowance. Old active ledger/guards/workers remain untouched. Freeze continuation clocks explicitly, retaining the original96h boundary rather than silently resetting it; new fullbody workers retain2h limits. The single transfer-family lock, deferred signals, owned cleanup/quarantine, original floors/RSS/resource meter and allowedonecurl-plusoneheavyworker ordering remain intact.

Exact all13 final proof cardinality, both identical primary receipts with frozen intended hashes, Terminal2 private PENDING/seal and current body/receipt/HTTP/resource/signal gates before and after persistence are required. Any new failure remains preserved; no automatic retry is admitted. Root will assign a different author to prepare the actual new plan/controller/common meter, followed by separate independent exact-version review and root admission. This diagnosis is not that future prelaunch PASS, source clearance, pipeline reproduction, replication or scientific completion.
"""
with MD.open('x') as stream:stream.write(text)
files.update(r['fixture_regular_sha256'])
for q in [CHECKER,RECEIPT,MD,REPORT,SPEC,Path(__file__)]:files[str(q)]=sha(q)
assert all(sha(q)==v for q,v in files.items())
save(SEAL,dict(schema='independent_validation13_resource_race_diagnosis_seal_v1',
    verdict=report['verdict'],file_sha256=files,regular_file_count=len(files),
    root30bodyprefix_receipt_evidence_bound_without_reference_traversal=True,
    no11body_reread_or_download_or_oldsuite=True,
    no_corrected_executor_prelaunch_or_continuation_execution_admission=True))
assert all(sha(q)==v for q,v in files.items())
print(json.dumps(dict(verdict=report['verdict'],md=str(MD),md_sha256=sha(MD),
    report=str(REPORT),report_sha256=sha(REPORT),spec=str(SPEC),spec_sha256=sha(SPEC),
    seal=str(SEAL),seal_sha256=sha(SEAL),file_count=len(files))))
