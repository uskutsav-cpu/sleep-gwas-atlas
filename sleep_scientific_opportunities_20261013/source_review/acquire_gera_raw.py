"""One reviewed source-only curl acquisition. No preprocessing or estimator."""
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

ROOT = Path(__file__).resolve().parent
PLAN = ROOT / 'GERA_RAW_ACQUISITION_PLAN.json'

def digest(path, deadline, internal_floor):
    md5, sha = hashlib.md5(), hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024**2), b''):
            if time.monotonic() > deadline:
                raise RuntimeError('VERIFICATION_WALL_TIME_GUARD')
            if shutil.disk_usage(ROOT).free < internal_floor:
                raise RuntimeError('VERIFICATION_INTERNAL_FLOOR')
            md5.update(block)
            sha.update(block)
    return md5.hexdigest(), sha.hexdigest()

def main():
    plan_bytes = PLAN.read_bytes()
    plan = json.loads(plan_bytes)
    review = json.loads((ROOT / 'GERA_RAW_ACQUISITION_REVIEW.json').read_text())
    plan_sha = hashlib.sha256(plan_bytes).hexdigest()
    script_sha = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    if (review['reviewed_plan_sha256'] != plan_sha or review['acquisition_admitted'] is not True
            or review['reviewed_script_sha256'] != script_sha):
        raise RuntimeError('FOCUSED_PLAN_ADMISSION_REQUIRED')
    out, part = Path(plan['destination']), Path(plan['checkpoint'])
    receipt_path = ROOT / 'GERA_RAW_ACQUISITION_RECEIPT.json'
    if out.exists() or part.exists() or receipt_path.exists():
        raise RuntimeError('EXISTING_SOURCE_CHECKPOINT_OR_RECEIPT_PRESERVED')
    resource, source = plan['resource'], plan['source']
    if shutil.disk_usage('/Volumes/Extreme SSD').free < resource['ssd_minimum_before_bytes']:
        raise RuntimeError('SSD_CAPACITY_GUARD')
    if shutil.disk_usage(ROOT).free < resource['internal_minimum_before_bytes']:
        raise RuntimeError('ACQUISITION_INTERNAL_GUARD')
    out.parent.mkdir(exist_ok=True)
    headers = out.parent / 'HTTP_headers_private.txt'
    stats = out.parent / 'curl_status_private.json'
    errors = out.parent / 'curl_stderr_private.txt'
    if any(p.exists() for p in [headers, stats, errors]):
        raise RuntimeError('PRIOR_TRANSPORT_EVIDENCE_PRESERVED')
    # curl writes only to an exclusively reserved checkpoint on the SSD.
    with part.open('xb'):
        pass
    command = ['/usr/bin/curl', '--disable', '--fail', '--silent', '--show-error', '--proto', '=https',
               '--connect-timeout', '30', '--max-time', str(resource['maximum_seconds']),
               '--max-filesize', str(resource['maximum_transfer_bytes']),
               '--header', 'Accept-Encoding: identity',
               '--header', 'If-Match: "'+source['etag']+'"',
               '--dump-header', str(headers), '--output', str(part),
               '--write-out', '%{json}', source['url']]
    stop, peak, started = None, 0, time.monotonic()
    deadline = started + resource['maximum_seconds']
    with stats.open('xb') as sout, errors.open('xb') as serr:
        process = subprocess.Popen(command, stdout=sout, stderr=serr)
        try:
            while process.poll() is None:
                rss = subprocess.run(['ps', '-p', str(process.pid), '-o', 'rss='],
                                     capture_output=True, text=True).stdout.strip()
                if rss:
                    peak = max(peak, int(rss)*1024)
                if peak > resource['curl_max_rss_bytes']:
                    stop = 'CURL_RSS_GUARD'
                elif shutil.disk_usage(ROOT).free < resource['internal_emergency_floor_bytes']:
                    stop = 'ACQUISITION_INTERNAL_FLOOR'
                elif part.stat().st_size > resource['maximum_retained_source_bytes']:
                    stop = 'RETAINED_SIZE_GUARD'
                elif time.monotonic() > deadline:
                    stop = 'WALL_TIME_GUARD'
                if stop:
                    break
                time.sleep(2)
        except Exception as error:
            stop = 'MONITOR_EXCEPTION_'+type(error).__name__
        finally:
            if process.poll() is None:
                process.terminate()
                try:
                    process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    process.kill()
                    process.wait()
    try:
        transport = json.loads(stats.read_text())
    except (ValueError, OSError):
        transport = {}
    actual_bytes = part.stat().st_size
    try:
        md5, sha = digest(part, deadline, resource['internal_emergency_floor_bytes'])
    except RuntimeError as error:
        md5, sha = None, None
        stop = str(error)
    raw_headers = headers.read_text() if headers.exists() else ''
    header_lines = raw_headers.splitlines()
    received_etags = [x.split(':', 1)[1].strip().strip('"') for x in header_lines if x.lower().startswith('etag:')]
    lengths = [x.split(':', 1)[1].strip() for x in header_lines if x.lower().startswith('content-length:')]
    verified = (process.returncode == 0 and stop is None and transport.get('http_code') == 200
                and actual_bytes == source['bytes'] and md5 == source['md5_from_yaml']
                and source['etag'] in received_etags and str(source['bytes']) in lengths
                and PLAN.read_bytes() == plan_bytes)
    if verified:
        with part.open('r+b') as stream:
            os.fsync(stream.fileno())
        if time.monotonic() > deadline:
            verified, stop = False, 'PUBLICATION_WALL_TIME_GUARD'
        elif shutil.disk_usage(ROOT).free < resource['internal_emergency_floor_bytes']:
            verified, stop = False, 'PUBLICATION_INTERNAL_FLOOR'
        else:
            # macOS BSD mv -n leaves any existing destination untouched.
            publish = subprocess.run(['mv', '-n', str(part), str(out)], capture_output=True, text=True)
            if publish.returncode != 0 or part.exists() or not out.exists():
                verified, stop = False, 'NON_CLOBBERING_PUBLICATION_NOT_COMPLETED'
    receipt = {'completed_utc': datetime.now(timezone.utc).isoformat(),
               'status': 'EXACT_SOURCE_ACQUISITION_PASS' if verified else 'FAILED_CHECKPOINT_PRESERVED',
               'plan_sha256': plan_sha, 'script_sha256': script_sha, 'url': source['url'],
               'path': str(out if verified else part), 'expected_bytes': source['bytes'],
               'actual_bytes': actual_bytes, 'expected_md5': source['md5_from_yaml'],
               'actual_md5': md5, 'actual_sha256': sha, 'exit_code': process.returncode,
               'http_code': transport.get('http_code'), 'received_etags': received_etags,
               'received_content_lengths': lengths, 'stop_reason': stop,
               'elapsed_seconds': time.monotonic()-started, 'curl_peak_observed_rss_bytes': peak,
               'internal_free_after_bytes': shutil.disk_usage(ROOT).free,
               'native_3GiB_guard_unchanged': True, 'source_only_acquisition': True,
               'original_sources_modified': False, 'h2_or_rg_runs': 0,
               'preprocessing_admitted': False, 'downstream_source_QC_status': 'NOT_TESTED'}
    with receipt_path.open('x') as stream:
        json.dump(receipt, stream, indent=2)
        stream.write('\n')
    print(json.dumps(receipt), flush=True)
    if not verified:
        raise RuntimeError('ACQUISITION_NOT_VERIFIED_CHECKPOINT_PRESERVED')

if __name__ == '__main__':
    main()
