#!/usr/bin/env python3
"""Record a complete current, read-only decoder B environment identity.

This is an identity census, not historical runtime recovery, decoder
qualification, source-body admission, or permission to modify the environment.
"""
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import resource
import stat
import subprocess
import time


PACKAGE = Path(__file__).resolve().parents[1]
ENVIRONMENT = Path('/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/FAILED-PARTIAL-2026-09-01-sleep_gwas_atlas/.ldsc-env')
PROFILE = PACKAGE / 'source_provenance/signed_reference_decoderB_runtime_identity_v1.json'
RECEIPT = PACKAGE / 'logs/signed_reference_decoderB_runtime_identity_audit_receipt_v1.json'
PYTHON = ENVIRONMENT / 'bin/python'
BINARY = ENVIRONMENT / 'bin/python3.9'
EXPECTED_BINARY_SHA = '7dffb088cd3027e48f0127ced6f206e06111abed3b9457e580b6ebbf593c10ba'
QUERY = "import sys,json,numpy,pandas,bitarray;print(json.dumps(dict(python=sys.version,executable=sys.executable,prefix=sys.prefix,base_prefix=sys.base_prefix,numpy=numpy.__version__,pandas=pandas.__version__,bitarray=bitarray.__version__)))"


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b''):
            digest.update(chunk)
    return digest.hexdigest()


def check_budget(start):
    if time.monotonic() - start > 300:
        raise RuntimeError('RUNTIME_IDENTITY_AUDIT_300_SECOND_LIMIT')
    if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss > 256 << 20:
        raise RuntimeError('RUNTIME_IDENTITY_AUDITOR_256_MIB_LIMIT')


def census(start):
    if ENVIRONMENT.is_symlink() or not ENVIRONMENT.is_dir():
        raise RuntimeError('REGULAR_READ_ONLY_ENVIRONMENT_REQUIRED')
    regular = {}
    links = {}
    directories = []
    total = 0
    for folder, children, files in os.walk(ENVIRONMENT, followlinks=False):
        directories.append(str(Path(folder)))
        for name in sorted(children + files):
            path = Path(folder) / name
            info = path.lstat()
            if stat.S_ISLNK(info.st_mode):
                resolved = path.resolve(strict=True)
                if not resolved.is_file() or resolved.is_symlink():
                    raise RuntimeError('REGULAR_RESOLVED_LINK_TARGET_REQUIRED: ' + str(path))
                links[str(path)] = dict(literal_link=os.readlink(path),
                    resolved_target=str(resolved), target_is_file=True,
                    resolved_sha256=sha(resolved))
            elif stat.S_ISREG(info.st_mode):
                if any(parent.is_symlink() for parent in path.parents):
                    raise RuntimeError('REGULAR_FILE_PARENT_REQUIRED: ' + str(path))
                regular[str(path)] = dict(bytes=info.st_size, sha256=sha(path))
                total += info.st_size
                if total > 1 << 30:
                    raise RuntimeError('READ_ONLY_RUNTIME_ONE_GIB_INVENTORY_BOUND')
            elif not stat.S_ISDIR(info.st_mode):
                raise RuntimeError('UNSUPPORTED_ENVIRONMENT_NODE: ' + str(path))
            check_budget(start)
    return regular, links, sorted(directories), total


def versions():
    result = subprocess.run([str(PYTHON), '-B', '-c', QUERY],
        env={**os.environ, 'PYTHONDONTWRITEBYTECODE': '1',
             'OPENBLAS_NUM_THREADS': '1', 'OMP_NUM_THREADS': '1',
             'MKL_NUM_THREADS': '1', 'NUMEXPR_NUM_THREADS': '1'},
        capture_output=True, text=True, check=True, timeout=60)
    value = json.loads(result.stdout)
    if not value['python'].startswith('3.9.23') or (
            value['numpy'], value['pandas'], value['bitarray']) != ('1.21.5', '1.3.3', '2.8.3'):
        raise RuntimeError('EXACT_DECODER_B_DECLARED_VERSIONS_REQUIRED')
    if value['executable'] != str(PYTHON) or value['prefix'] != str(ENVIRONMENT):
        raise RuntimeError('EXACT_DECODER_B_LOGICAL_ENVIRONMENT_REQUIRED')
    return value


def write_new(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write('\n')
        stream.flush()
        os.fsync(stream.fileno())


def main():
    if PROFILE.exists() or RECEIPT.exists():
        raise RuntimeError('EXISTING_RUNTIME_IDENTITY_EVIDENCE_PRESERVED_NO_OVERWRITE')
    if os.readlink(PYTHON) != 'python3.9' or BINARY.is_symlink() or sha(BINARY) != EXPECTED_BINARY_SHA:
        raise RuntimeError('EXACT_INHERITED_B_BINARY_AND_LITERAL_LINK_REQUIRED')
    start = time.monotonic()
    before = census(start)
    first = versions()
    last = versions()
    after = census(start)
    if before != after or first != last:
        raise RuntimeError('COMPLETE_ENVIRONMENT_OR_VERSIONS_CHANGED_DURING_READ_ONLY_AUDIT')
    check_budget(start)
    profile = dict(status='READ_ONLY_RUNTIME_CANDIDATE_MATCHES_DECLARED_DECODER_B_VERSIONS',
        created_utc=datetime.now(timezone.utc).isoformat(), environment_root=str(ENVIRONMENT),
        command=[str(PYTHON), '-B', '-c', QUERY], versions=first,
        versions_before=first, versions_after=last, regular_files=before[0],
        symlinks=before[1], directories=before[2], regular_file_count=len(before[0]),
        symlink_count=len(before[1]), regular_file_bytes=before[3],
        all_before_after_identities_unchanged=True, complete_environment_census=True,
        inherited_resolved_python_binary_sha256=EXPECTED_BINARY_SHA,
        pycache_writes_disabled=True, BLAS_threads=1,
        other_project_files_changed=False, data_body_reads=0, worker_fit_calls=0,
        runtime_installs_or_copies=0, historical_per_trait_runtime_attestation_recovered=False,
        exact_historical_binary_recovery_claimed=False,
        qualification='Complete current environment census and declared versions only. System libraries, historical execution identity and scientific decoder/source qualification remain separate. Read-only use requires fresh current inventory verification and independent review.',
        elapsed_seconds=time.monotonic()-start,
        peak_auditor_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    write_new(PROFILE, profile)
    receipt = dict(schema='root_decoderB_complete_current_read_only_runtime_identity_audit_v1',
        recorded_utc=datetime.now(timezone.utc).isoformat(), producer_script=str(Path(__file__)),
        producer_script_sha256=sha(__file__), profile=str(PROFILE), profile_sha256=sha(PROFILE),
        regular_file_count=len(before[0]), symlink_count=len(before[1]),
        regular_file_bytes=before[3], complete_two_censuses_equal=True,
        versions_before_after_equal=True, source_or_genotype_body_access=False,
        real_decoder_execution_or_estimator_calls=False,
        independent_runtime_profile_review_pass=False,
        source_stage_execution_admitted=False,
        elapsed_seconds=time.monotonic()-start)
    write_new(RECEIPT, receipt)
    print(json.dumps(dict(status=profile['status'], profile_sha256=sha(PROFILE),
        receipt_sha256=sha(RECEIPT), regular_files=len(before[0]),
        symlinks=len(before[1]), bytes=before[3], elapsed_seconds=receipt['elapsed_seconds'])))


if __name__ == '__main__':
    main()
