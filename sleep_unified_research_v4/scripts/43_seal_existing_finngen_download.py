#!/usr/bin/env python3
"""Explicit audit of a completed transfer whose final path check failed; no retry."""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil

ROOT = Path(__file__).resolve().parents[2]
PACKAGE = ROOT / 'sleep_unified_research_v4'
SSD = Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4')


def sha(path, algorithm='sha256'):
    h = hashlib.new(algorithm)
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()


def main():
    before = PACKAGE / 'logs/finngen_R13_insomnia_acquisition_receipt_v4.json'
    previous = json.loads(before.read_text())
    assert previous['status'] == 'FAILED_PRESERVED_NO_RETRY'
    assert previous['returncode'] == 0
    assert previous['stop_reason'] == "AttributeError: 'str' object has no attribute 'open'"
    assert previous['teardown']['remaining_group_members'] == []
    folder = SSD / 'new_source_feasibility/finngen_R13_F5_INSOMNIA'
    partial = folder / 'finngen_R13_F5_INSOMNIA.gz.partial'
    final = folder / 'finngen_R13_F5_INSOMNIA.gz'
    receipt_path = PACKAGE / 'logs/finngen_R13_insomnia_acquisition_verified_v4_2.json'
    assert partial.is_file() and not final.exists() and not receipt_path.exists()
    assert partial.stat().st_size == previous['expected_size'] == 809346932
    assert sha(partial, 'md5') == previous['expected_md5'] == 'f16c21acbf8b9ebfeb0f26a19f0ecfe3'
    header = folder / 'download_headers.txt'
    assert sha(header) == previous['headers_sha256']
    assert 'x-goog-generation: 1777989563097164' in header.read_text().lower()
    assert all(sha(p) == h for p, h in previous['bound_sources'].items())
    internal, ssd = shutil.disk_usage('/System/Volumes/Data').free, shutil.disk_usage(SSD).free
    assert internal >= 3 * 1024**3 and ssd >= 5 * 1024**3
    body_sha = sha(partial)
    evidence = {'recorded_utc': datetime.now(timezone.utc).isoformat(),
                'scope': 'EXPLICIT_REVIEW_OF_EXISTING_COMPLETE_TRANSFER_NO_REDOWNLOAD',
                'original_failed_receipt': str(before), 'original_failed_receipt_sha256': sha(before),
                'original_failure_preserved': True, 'new_verification_code_sha256': sha(Path(__file__)),
                'source_bound_hashes_unchanged': True, 'headers_sha256': sha(header),
                'generation': '1777989563097164', 'actual_size': partial.stat().st_size,
                'actual_md5': sha(partial, 'md5'), 'actual_sha256': body_sha,
                'source_path': str(final), 'status': 'EXACT_GENERATION_SIZE_MD5_ACQUIRED_SEALED_AFTER_EXPLICIT_AUDIT',
                'gzip_full_stream_verified': False, 'scientific_source_admitted': False,
                'new_network_requests': 0, 'internal_free_bytes': internal, 'ssd_free_bytes': ssd}
    partial.rename(final)
    assert sha(final) == body_sha
    receipt_path.write_text(json.dumps(evidence, indent=2) + '\n')
    (folder / receipt_path.name).write_bytes(receipt_path.read_bytes())
    print(json.dumps({k: evidence[k] for k in ('status', 'actual_size', 'actual_sha256', 'new_network_requests')}))


if __name__ == '__main__':
    main()
