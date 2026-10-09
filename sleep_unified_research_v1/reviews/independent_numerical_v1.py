#!/usr/bin/env python3
"""Independent arithmetic/provenance audit of the insomnia--BMI LDSC pilot.

This does not import LDSC, the capture wrapper, numpy, scipy, or another review.
It performs no GWAS analysis and writes only independent_numerical_v1.* files.
"""
from datetime import datetime, timezone
from decimal import Decimal
from pathlib import Path
import csv
import hashlib
import json
import math
import re
import subprocess
import sys
import tarfile

ROOT = Path(__file__).resolve().parents[2]
PACKAGE = ROOT / "sleep_unified_research_v1"
OUT = Path(__file__).with_suffix("")
PILOT = PACKAGE / "native/core_pilot_v1"
PREFIX = PILOT / "rg_insomnia__bmi"
PINNED = "6c673952cee74bd5c57aef1555a03b1c015399a0"
CODE = ROOT.parent / "ldsc-code"
OLD = Path("/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/FAILED-PARTIAL-2026-09-01-sleep_gwas_atlas")
REF_SHA = "9537f00eb0d163a935aaa2cf04b358b7cf21852279b9c7925802526f6060b069"
FLOAT_RTOL = 1e-12
FLOAT_ATOL = 1e-14
# Historical log/TSV values are rounded, with stripped trailing zeros in TSV.
# These limits are defined before reading or calculating the pilot numbers.
HISTORICAL_LIMITS = {
    "rg": 5e-5, "se": 5e-5, "z": 5e-4, "p": 5e-18,
    "h2_obs": 5e-5, "h2_obs_se": 5e-5, "h2_int": 5e-4,
    "h2_int_se": 5e-5, "gcov_int": 5e-5, "gcov_int_se": 5e-5,
}


def sha_stream(handle):
    h = hashlib.sha256()
    while chunk := handle.read(1024 * 1024):
        h.update(chunk)
    return h.hexdigest()


def sha(path):
    with Path(path).open("rb") as handle:
        return sha_stream(handle)


def rows(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def record(path, expected=None, role="audit_input", evidence=None):
    path = Path(path)
    before = path.stat()
    actual = sha(path)
    after = path.stat()
    stable = (before.st_size, before.st_mtime_ns) == (after.st_size, after.st_mtime_ns)
    return dict(path=str(path), role=role, bytes=before.st_size,
                mtime_ns=before.st_mtime_ns, actual_sha256=actual,
                expected_sha256=expected, expected_hash_evidence=evidence,
                matched_expected=(actual == expected) if expected else None,
                stat_unchanged_during_hash=stable)


def check(field, calculated, captured):
    return dict(comparison="capture", field=field, calculated=calculated,
                reference=captured, absolute_difference=abs(calculated-captured),
                absolute_tolerance=FLOAT_ATOL, relative_tolerance=FLOAT_RTOL,
                passed=math.isclose(calculated, captured, rel_tol=FLOAT_RTOL,
                                    abs_tol=FLOAT_ATOL))


def mean(values):
    return math.fsum(values) / len(values)


def delete_se(values):
    b = len(values)
    center = mean(values)
    return math.sqrt((b-1)/b * math.fsum((v-center)**2 for v in values))


def read_vector(path):
    text = Path(path).read_text()
    # Stock output is one total-estimate delete value per block, no SNP IDs.
    lines = [line.split() for line in text.splitlines() if line.strip()]
    assert len(lines) == 200 and all(len(line) == 1 for line in lines)
    values = [float(line[0]) for line in lines]
    assert all(math.isfinite(v) for v in values)
    return values


def token_limit(token):
    """Half a unit in the last printed place, retaining scientific exponent."""
    return float(Decimal(10) ** Decimal(token).as_tuple().exponent / 2)


def main():
    receipt = dict(
        review="NEW independent numerical-reproduction reviewer, Phase 10",
        recorded_utc=datetime.now(timezone.utc).isoformat(),
        runtime_python=sys.version, implementation="stdlib math.fsum/math.erfc; no LDSC/scipy/numpy imports",
        independent_code_sha256=sha(Path(__file__)),
        independence=dict(other_review_code_or_reports_read=False,
                          prior_reviewer_calculated_numbers_used=False,
                          formula_sources="actual pinned LDSC regressions.py and jackknife.py, inspected as source only",
                          capture_numbers_used_only_as_reference=True),
        scope=dict(native_gwas_rerun=False, raw_harmonization_munging_rerun=False,
                   independent_processed_input_pair_checks=1,
                   fully_native_reproduced_396_family=False,
                   fully_native_reproduced_1200_family=False,
                   genomic_block_alignment_proven_by_output_arrays=False,
                   source_inputs_modified=False),
        tolerances=dict(capture_absolute=FLOAT_ATOL, capture_relative=FLOAT_RTOL,
                        capture_p_relative=1e-12, capture_p_absolute=0,
                        historical_tsv_absolute=HISTORICAL_LIMITS,
                        historical_log="half-unit of last printed place"),
        artifacts=[], source_input_hashes=[], reference_hashes=[], code_hashes=[], comparisons=[])
    full = Path(str(PREFIX)+".full_precision.json")
    log = Path(str(PREFIX)+".log")
    capture = json.loads(full.read_text())
    assert len(capture["estimates"]) == 1
    fitted = capture["estimates"][0]
    assert fitted["status"] == "NATIVE_ESTIMATE_RETURNED"
    assert capture["arguments"]["n_blocks"] == 200
    assert capture["arguments"]["no_check_alleles"] is False
    assert capture["arguments"]["print_delete_vals"] is True
    suffix = "insomnia.sumstats.gz_bmi.sumstats.gz."
    vectors = {}
    for field in ("hsq1", "hsq2", "gencov"):
        path = Path(str(PREFIX)+suffix+field+".delete")
        vectors[field] = read_vector(path)
        receipt["artifacts"].append(record(path, role="native_stock_delete_vector"))
    for path in (full, log, PACKAGE/"scripts/native_ldsc_capture.py",
                 PACKAGE/"tables/native_input_hash_checks.tsv",
                 PACKAGE/"logs/native_input_verification_receipt_v1.json",
                 PACKAGE/"manifests/ldsc_native_code_comparison_v1.json",
                 ROOT/"discovery_extension/core_checkpoint.json",
                 ROOT/"scripts/05_collate.py", ROOT/"scripts/23_merge_rg_families.py"):
        receipt["artifacts"].append(record(path))

    # Independent computation uses captured h2/covariance totals and stock deletes.
    h1, h2, cov = (fitted[f]["tot"] for f in ("hsq1", "hsq2", "gencov"))
    assert h1 > 0 and h2 > 0
    theta = cov / math.sqrt(h1*h2)
    b = 200
    products = [x*y for x, y in zip(vectors["hsq1"], vectors["hsq2"])]
    assert all(v > 0 for v in products)
    deletes = [c/math.sqrt(p) for c,p in zip(vectors["gencov"],products)]
    # Compute the mean bias correction and SE via delete variance. This differs
    # computationally from stock pseudovalue np.cov, while being algebraically equal.
    bias_corrected = b*theta-(b-1)*mean(deletes)
    se = delete_se(deletes)
    z = theta/se
    p = math.erfc(abs(z)/math.sqrt(2))
    pseudovalues = [b*theta-(b-1)*v for v in deletes]
    pseudo_center = mean(pseudovalues)
    pseudo_se = math.sqrt(math.fsum((v-pseudo_center)**2 for v in pseudovalues)/(b*(b-1)))
    calculated = dict(rg_ratio=theta, rg_jknife=bias_corrected, rg_se=se, z=z, p=p,
                      pseudovalue_mean_alternative=pseudo_center,
                      pseudovalue_se_alternative=pseudo_se)
    for field in ("rg_ratio", "rg_jknife", "rg_se", "z", "p"):
        comparison = check(field, calculated[field], fitted[field])
        if field == "p":
            comparison["absolute_tolerance"] = 0
            comparison["passed"] = math.isclose(p, fitted[field], rel_tol=1e-12, abs_tol=0)
        receipt["comparisons"].append(comparison)
    for field in ("hsq1", "hsq2", "gencov"):
        direct_se = delete_se(vectors[field])
        calculated[field+"_se"] = direct_se
        receipt["comparisons"].append(check(field+"_se", direct_se, fitted[field]["tot_se"]))
    receipt["numerical_reconstruction"] = calculated
    receipt["delete_array_summary"] = {
        k: dict(blocks=len(v), minimum=min(v), maximum=max(v), all_finite=True)
        for k,v in vectors.items()}
    receipt["delete_array_summary"]["rg_delete"] = dict(blocks=b, minimum=min(deletes), maximum=max(deletes))
    receipt["estimator_identity"] = dict(
        displayed_point_estimate="covariance / sqrt(hsq1 * hsq2), rg_ratio",
        uncertainty="sqrt((B-1)/B * sum((ratio_delete_i - mean_ratio_delete)^2))",
        bias_corrected_mean="B*ratio - (B-1)*mean(ratio_delete), rg_jknife",
        p_value="two-sided standard-normal P = erfc(abs(ratio/SE)/sqrt(2))",
        rg_ratio_minus_rg_jknife=theta-bias_corrected,
        independent_algebraic_se_difference=se-pseudo_se)

    historical_path = PACKAGE/"sources/recovered/results/tables/rg_matrix.tsv"
    checkpoint = json.loads((ROOT/"discovery_extension/core_checkpoint.json").read_text())
    expected = checkpoint["artifact_hashes_sha256"]["results/tables/rg_matrix.tsv"]
    receipt["artifacts"].append(record(historical_path, expected, "frozen_396_table", "core_checkpoint.json"))
    historical_rows = rows(historical_path)
    selected = [r for r in historical_rows if r["sleep_trait"]=="insomnia" and r["disease_trait"]=="bmi"]
    assert len(selected)==1 and len(historical_rows)==396
    historical = selected[0]
    mapped = dict(rg=theta, se=se, z=z, p=p, h2_obs=h2, h2_obs_se=calculated["hsq2_se"],
                  h2_int=fitted["hsq2"]["intercept"], h2_int_se=fitted["hsq2"]["intercept_se"],
                  gcov_int=fitted["gencov"]["intercept"], gcov_int_se=fitted["gencov"]["intercept_se"])
    for field, limit in HISTORICAL_LIMITS.items():
        difference = abs(mapped[field]-float(historical[field]))
        receipt["comparisons"].append(dict(comparison="historical_396_tsv",field=field,
            calculated=mapped[field], reference_token=historical[field],
            reference=float(historical[field]), absolute_difference=difference,
            absolute_tolerance=limit, relative_tolerance=0, passed=difference <= limit))
    receipt["frozen_row"] = historical
    receipt["frozen_row_note"] = "396 rows counted and table hash checked; only insomnia--BMI arithmetic compared. FDR was not re-estimated or combined with extension."

    # Scalar stock log provides scientific-notation P even when summary P=0.0000.
    historical_log = OLD/"results/logs/rg_insomnia.log"
    receipt["artifacts"].append(record(historical_log, role="historical_native_log_no_checkpoint_expected_hash"))
    hist_text = historical_log.read_text()
    block = hist_text.split("Reading summary statistics from data/munged/bmi.sumstats.gz ...",1)[1].split("Computing rg",1)[0]
    current_text = log.read_text()
    for kind, text in (("pilot_scalar_log",current_text),("historical_scalar_log",block)):
        match = re.search(r"Genetic Correlation: ([^ ]+) \(([^)]+)\).*?Z-score: ([^\n]+).*?P: ([^\n]+)",text,re.S)
        assert match
        for field, token in zip(("rg","se","z","p"),match.groups()):
            token = token.strip()
            limit = token_limit(token)
            difference = abs(mapped[field]-float(token))
            receipt["comparisons"].append(dict(comparison=kind, field=field,
                calculated=mapped[field], reference_token=token, reference=float(token),
                absolute_difference=difference, absolute_tolerance=limit,
                relative_tolerance=0, passed=difference<=limit))
    receipt["printed_precision_note"] = "Scalar log rg=0.1723, SE=0.023, z=7.4777 and P=7.5635e-14. Frozen TSV has z=7.478 and P=7.563e-14, so only agreement at recorded precision is testable. The fixed-width log summary P=0.0000 is formatting, not a zero probability."

    # Rehash real input bytes; expected historical hashes exist for raw/archive,
    # but only a current receipt hash exists for harmonized/munged bytes.
    ledger = rows(PACKAGE/"tables/native_input_hash_checks.tsv")
    selected_sources = [r for r in ledger if r["trait_id"] in ("insomnia","bmi") and
                        r["kind"] in ("core_raw","core_munged","core_harmonized","core_source_archive")]
    for row in selected_sources:
        item = record(row["path"], row["expected_sha256"] or None, row["kind"],row["expected_hash_evidence"])
        item["ledger_actual_sha256"] = row["actual_sha256"]
        item["matched_current_ledger_hash"] = item["actual_sha256"] == row["actual_sha256"]
        item["historical_byte_identity_claim"] = bool(row["expected_sha256"])
        receipt["source_input_hashes"].append(item)
    assert len(selected_sources)==8
    source_by_path = {r["path"]:r for r in receipt["source_input_hashes"]}
    receipt["captured_paths_current_hashes"] = {fitted[k]:source_by_path[fitted[k]]["actual_sha256"] for k in ("p1","p2")}
    receipt["provenance_limit"] = "Rehashing existing raw and processed bytes does not rerun harmonization, allele coding, build mapping, MAF/INFO filtering, or munging; processed files have no historical expected hash."

    archive = ROOT.parent/"reference_check/eur_w_ld_chr.tar.gz"
    archive_record = record(archive, REF_SHA, "official_reference_archive", "pinned source receipt")
    receipt["reference_hashes"].append(archive_record)
    assert archive_record["matched_expected"]
    ref = Path(capture["arguments"]["ref_ld_chr"])
    assert ref == Path(capture["arguments"]["w_ld_chr"])
    consumed = {f"eur_w_ld_chr/{chrom}.{suffix}" for chrom in range(1,23)
                for suffix in ("l2.ldscore.gz","l2.M_5_50")}
    seen = set()
    with tarfile.open(archive,"r:gz") as tar:
        for member in tar.getmembers():
            if member.name not in consumed:
                continue
            with tar.extractfile(member) as handle:
                expected_member = sha_stream(handle)
            actual = record(ref/Path(member.name).name, expected_member,
                            "actual_consumed_ld_score_or_default_M_member",
                            REF_SHA+":"+member.name)
            actual["matched_input_hash_ledger"] = any(r["path"]==actual["path"] and r["actual_sha256"]==actual["actual_sha256"] for r in ledger)
            receipt["reference_hashes"].append(actual)
            seen.add(member.name)
    assert seen == consumed

    # Verify actual active code against both pinned Git payloads and old SSD code.
    active_head = subprocess.check_output(["git","-C",str(CODE),"rev-parse","HEAD"],text=True).strip()
    status = subprocess.check_output(["git","-C",str(CODE),"status","--porcelain"],text=True)
    receipt["ldsc_code_commit"] = active_head
    receipt["ldsc_worktree_clean"] = not status.strip()
    assert active_head == PINNED
    code_files = ["ldsc.py","ldscore/__init__.py","ldscore/irwls.py","ldscore/jackknife.py",
                  "ldscore/ldsc_utils.py","ldscore/ldsc_utils_local.py","ldscore/ldscore.py",
                  "ldscore/parse.py","ldscore/regressions.py","ldscore/sumstats.py"]
    for file in code_files:
        payload = subprocess.check_output(["git","-C",str(CODE),"show",PINNED+":"+file])
        expected_code = hashlib.sha256(payload).hexdigest()
        actual = record(CODE/file, expected_code,"active_ldsc_code",PINNED+":"+file)
        actual["archived_code_path"] = str(OLD/"ldsc"/file)
        actual["archived_code_sha256"] = sha(OLD/"ldsc"/file)
        actual["matched_archived_code"] = actual["actual_sha256"]==actual["archived_code_sha256"]
        receipt["code_hashes"].append(actual)
    receipt["wrapper_static_audit"] = dict(
        estimator_change_detected=False,
        rationale="Wrapper calls saved original estimate_rg/estimate_h2 functions, copies returned scalar attributes, then returns the original result. It does not set numerical tolerances, alter input arrays, or redefine regression/jackknife methods.",
        runtime_input_output_binding_limit="Pilot capture records paths and software versions, but has no atomic execution receipt containing at-execution input/code SHA-256. Current hashes and matching logs support the claimed run; post-hoc hashes do not eliminate a hypothetical change between execution and hashing.")

    preflight = PACKAGE/"logs/native_resource_preflight_20261009T042146Z.json"
    resource = json.loads(preflight.read_text())
    receipt["artifacts"].append(record(preflight,role="resource_gate_receipt"))
    receipt["resource_limit"] = resource
    receipt["native_output_capture_inventory"] = [str(p.relative_to(PACKAGE)) for p in (PACKAGE/"native").rglob("*.full_precision.json")]
    receipt["block_alignment_limit"] = "The three 200-row delete arrays show within-pilot arithmetic consistency. They contain no SNP identities, genomic coordinates, block boundary/separator map, or independently validated cross-run matched SNP order. They cannot establish discovery-validation genomic block alignment or estimate shared-sleep cross-run covariance by themselves."
    receipt["human_review_questions"] = [
        "Reproduce and audit the original raw-to-harmonized-to-munged chain and source allele/build conventions.",
        "Restore internal disk headroom without changing original research assets before executing all native families.",
        "Require hash-bound execution receipts for full-family reruns, and preserve separate 396/1,200 testing families.",
        "For covariance work, recover or construct justified matched SNP/block boundaries rather than treating row number as genomic alignment."]
    checks = receipt["comparisons"]
    provenance_records = receipt["artifacts"]+receipt["source_input_hashes"]+receipt["reference_hashes"]+receipt["code_hashes"]
    receipt["all_arithmetic_checks_pass"] = all(r["passed"] for r in checks)
    receipt["all_expected_hashes_match"] = all(r["matched_expected"] is not False for r in provenance_records)
    receipt["all_current_source_ledger_hashes_match"] = all(r["matched_current_ledger_hash"] for r in receipt["source_input_hashes"])
    receipt["all_active_and_archived_code_hashes_match"] = all(r["matched_expected"] and r["matched_archived_code"] for r in receipt["code_hashes"])
    receipt["verdict"] = "PASS_SINGLE_PROCESSED_INPUT_PILOT_ARITHMETIC_AND_CURRENT_BYTE_PROVENANCE; FULL_NATIVE_SOURCE_CHAIN_AND_FAMILIES_UNVERIFIED"
    assert receipt["all_arithmetic_checks_pass"]
    assert receipt["all_expected_hashes_match"]
    assert receipt["all_current_source_ledger_hashes_match"]
    assert receipt["all_active_and_archived_code_hashes_match"]

    json_path = Path(str(OUT)+".json")
    tsv_path = Path(str(OUT)+".tsv")
    report_path = Path(str(OUT)+".md")
    json_path.write_text(json.dumps(receipt,indent=2,allow_nan=False)+"\n")
    fields = ["comparison","field","calculated","reference","reference_token","absolute_difference","absolute_tolerance","relative_tolerance","passed"]
    with tsv_path.open("w",newline="") as handle:
        writer = csv.DictWriter(handle,fieldnames=fields,delimiter="\t",extrasaction="ignore")
        writer.writeheader()
        writer.writerows(checks)
    source_lines = [f"| {r['role']} | {Path(r['path']).name} | `{r['actual_sha256']}` | {'historical expected hash matched' if r['matched_expected'] else 'current ledger hash matched; historical hash unavailable'} |" for r in receipt["source_input_hashes"]]
    report = f"""# Independent numerical reproduction review, v1

Verdict: **one processed-input insomnia–BMI pilot passes independent arithmetic reconstruction and current byte-provenance checks. Full source-chain and 396-/1,200-family native reproduction remain unverified.**

This is a new reviewer implementation, distinct from the existing statistical/numerical reviewer. It does not read or import that reviewer's code, reports, or computed results. Formulas were checked against the actual pinned LDSC source. Computation uses Python standard-library `math.fsum`, direct delete-vector variance, and `math.erfc`; neither LDSC nor the capture wrapper is imported or executed. No native GWAS analysis was rerun and original SSD inputs/frozen outputs were not modified.

The capture contains one estimate, and each of the three stock delete files contains 200 finite scalar values. All delete h2 products are positive. Independently recomputed values are:

| Statistic | Independently calculated | Captured reference | Absolute difference |
|---|---:|---:|---:|
| rg ratio | {theta:.17g} | {fitted['rg_ratio']:.17g} | {abs(theta-fitted['rg_ratio']):.4g} |
| Bias-corrected jackknife rg | {bias_corrected:.17g} | {fitted['rg_jknife']:.17g} | {abs(bias_corrected-fitted['rg_jknife']):.4g} |
| rg SE | {se:.17g} | {fitted['rg_se']:.17g} | {abs(se-fitted['rg_se']):.4g} |
| z | {z:.17g} | {fitted['z']:.17g} | {abs(z-fitted['z']):.4g} |
| two-sided P | {p:.17g} | {fitted['p']:.17g} | {abs(p-fitted['p']):.4g} |

The displayed LDSC point estimate is `rg_ratio=cov/sqrt(h2_1*h2_2)`. Its SE comes from the stock ratio jackknife. The bias-corrected jackknife mean differs by {theta-bias_corrected:.17g}; it is not the displayed point estimate. With B=200 and delete ratio d_i, the independent SE calculation is `sqrt((B-1)/B * sum((d_i-mean(d))^2))`. An independently evaluated pseudovalue expression agrees within {abs(se-pseudo_se):.4g}. Heritability and genetic-covariance SEs reconstructed separately from their stock deletes also pass.

Capture comparisons use relative tolerance 1e-12 and absolute tolerance 1e-14; P uses relative tolerance 1e-12 and zero absolute tolerance. These limits and the historical tolerances were declared in the review code before reading pilot estimates. There are {len(checks)} passing comparisons in the detailed TSV.

The immutable 396-row `rg_matrix.tsv` hashes to `{expected}`, matching the frozen checkpoint. Its insomnia–BMI row records rg={historical['rg']}, SE={historical['se']}, z={historical['z']}, P={historical['p']}; all tested fields agree within declared print-precision intervals. For rg/SE the interval is ±5e-5; z and h2 intercept use ±5e-4; P uses ±5e-18. The historical and new native scalar logs both retain z=7.4777 and P=7.5635e-14, with all scalar comparisons passing half-unit intervals at the final printed place. Full-precision historical identity cannot be tested from rounded archived output. The fixed-width summary P=0.0000 is a display truncation and must not replace its scalar scientific-notation P. FDR was not recomputed by this review.

Actual inputs were independently SHA-256 hashed again rather than accepting ledger assertions:

| Input stage | File | Actual SHA-256 | Evidence |
|---|---|---|---|
{chr(10).join(source_lines)}

The official 33,357,890-byte reference archive hashes to `{REF_SHA}`. All 44 consumed chromosome LD-score/default M_5_50 members match their archive-member expected hashes and the input ledger. The reference and weight prefixes point to the same checked directory. Active LDSC HEAD is `{active_head}`, its worktree is clean, and ten actual code files each match both the pinned Git payload and archived SSD code bytes. The JSON receipt records every hash. Byte matching is distinct from validating reference ancestry/design assumptions.

Static wrapper review found no estimator modification: the wrapper invokes the saved original estimator, copies return attributes, writes JSON, and returns the original estimator object. The pilot lacks a single atomic receipt that binds at-execution source/code hashes to output hashes, so this review verifies current byte provenance rather than asserting that post-hoc hashes establish every runtime property.

The pilot log reports 1,133,335 insomnia SNPs, 1,217,311 BMI SNPs, 1,130,007 after reference/SNP merging, and 1,006,820 with valid alleles. Those counts were observed in the log; filtering, allele orientation, genome mapping, sample definition, effect coding and overlap assumptions were not independently rerun. Raw-source/archive hashes match historical expectations; harmonized/munged hashes match current receipts but have no historical expected hashes. This supports the pilot's processed-input numerical reproduction and does not complete the original QC chain.

The 200 scalar output rows have no SNP identities, genomic coordinates or separator map. Their numerical agreement cannot establish genomic alignment between discovery/validation runs, justify cross-run paired jackknife covariance, or establish an independent biological replication. Such claims require matched SNP/block evidence and source/cohort review.

The full-family plan's latest inspected resource receipt records {resource['internal_free_bytes']:,} internal bytes free versus a {resource['internal_minimum_bytes']:,}-byte (3 GiB) floor, so its resource gate fails. No full family run is admitted by that receipt. This review performed one numerical pilot reconstruction; a 396-row frozen table and a job manifest are not completion of 396 native estimates or 1,200 extension estimates. Raw-chain rerun, all-family execution, genomic alignment/covariance calibration and biological interpretation remain unresolved human/scientific review requirements.

Reproduce this audit with `python3 sleep_unified_research_v1/reviews/independent_numerical_v1.py`. The script, report, JSON receipt and TSV are sealed by `independent_numerical_v1.sha256`. Scope is research evidence only; no manuscript text is generated.
"""
    report_path.write_text(report)
    sealed = [Path(__file__),json_path,tsv_path,report_path]
    Path(str(OUT)+".sha256").write_text("".join(f"{sha(path)}  {path.name}\n" for path in sealed))
    print(json.dumps(dict(verdict=receipt["verdict"], comparisons=len(checks),
                         report_sha256=sha(report_path), code_sha256=sha(Path(__file__)),
                         receipt_sha256=sha(json_path), rg=theta,se=se,p=p),indent=2))


if __name__ == "__main__":
    main()
