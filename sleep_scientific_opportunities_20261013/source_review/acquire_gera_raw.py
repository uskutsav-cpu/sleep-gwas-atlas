"""One reviewed source-only curl acquisition. No preprocessing or estimator."""
from datetime import datetime, timezone
import errno
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent
PLAN = ROOT / 'GERA_RAW_ACQUISITION_PLAN.json'

def persist_receipt(receipt, internal_path, ssd_path):
    """Persist terminal evidence on SSD before the optional internal mirror."""
    payload = json.dumps(receipt, indent=2) + '\n'
    with ssd_path.open('x') as stream:
        stream.write(payload)
        stream.flush()
        os.fsync(stream.fileno())
    try:
        with internal_path.open('x') as stream:
            stream.write(payload)
            stream.flush()
            os.fsync(stream.fileno())
    except OSError as error:
        if error.errno != errno.ENOSPC:
            raise
        return False
    return True

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
    ssd_receipt_path = part.parent / receipt_path.name
    if any(p.exists() for p in [out, part, receipt_path, ssd_receipt_path]):
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
    persist_receipt(receipt, receipt_path, ssd_receipt_path)
    print(json.dumps(receipt), flush=True)
    if not verified:
        raise RuntimeError('ACQUISITION_NOT_VERIFIED_CHECKPOINT_PRESERVED')

def range_header_check(raw, source, offset):
    """Return None until a complete final header block, otherwise exact verdict."""
    text = raw.replace('\r\n', '\n')
    blocks = text.split('\n\n')[:-1]
    final = []
    for block in blocks:
        if not block.startswith('HTTP/'):
            continue
        tokens = block.splitlines()[0].split()
        if len(tokens) >= 2 and len(tokens[1]) == 3 and tokens[1].isdigit() and tokens[1].startswith('1'):
            continue
        final.append(block)
    if not final:
        return None
    lines = final[-1].splitlines()
    fields = {}
    for line in lines[1:]:
        if ':' in line:
            key, value = line.split(':', 1)
            fields.setdefault(key.lower(), []).append(value.strip())
    expected = {'etag': ['"'+source['etag']+'"'],
                'content-length': [str(source['bytes']-offset)],
                'content-range': ['bytes '+str(offset)+'-'+str(source['bytes']-1)+'/'+str(source['bytes'])]}
    status = lines[0].split()
    valid = (len(status) >= 2 and status[1] == '206'
             and all(fields.get(key) == value for key, value in expected.items())
             and fields.get('content-encoding', ['identity']) == ['identity'])
    return {'valid': valid, 'status_line': lines[0], 'fields': fields}


def resume_main(plan_path):
    """Recover the owned prefix using a separate suffix, never overwrite it."""
    plan_path = Path(plan_path)
    plan_bytes = plan_path.read_bytes()
    plan = json.loads(plan_bytes)
    source, limits, saved = plan['source'], plan['resource'], plan['prefix']
    plan_sha = hashlib.sha256(plan_bytes).hexdigest()
    script_sha = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    review = json.loads((ROOT / plan['review']).read_text())
    if (review['reviewed_plan_sha256'] != plan_sha or review['reviewed_script_sha256'] != script_sha
            or review['resume_admitted'] is not True):
        raise RuntimeError('EXACT_RESUME_ADMISSION_REQUIRED')
    prefix, suffix = Path(saved['path']), Path(plan['suffix'])
    assembled, out = Path(plan['assembled']), Path(plan['destination'])
    receipt_path = ROOT / plan['receipt']
    ssd_receipt_path = suffix.parent / receipt_path.name
    if any(p.exists() for p in [suffix, assembled, out, receipt_path, ssd_receipt_path]):
        raise RuntimeError('PRIOR_PREFIX_SUFFIX_SOURCE_OR_RECEIPT_PRESERVED')
    started = time.monotonic()
    deadline = started + limits['maximum_seconds']
    floor = limits['internal_emergency_floor_bytes']
    if shutil.disk_usage(ROOT).free < limits['internal_minimum_before_bytes']:
        raise RuntimeError('RESUME_INTERNAL_GUARD')
    if shutil.disk_usage('/Volumes/Extreme SSD').free < limits['ssd_minimum_before_bytes']:
        raise RuntimeError('RESUME_SSD_GUARD')
    prefix_md5, prefix_sha = digest(prefix, deadline, floor)
    if prefix.stat().st_size != saved['bytes'] or prefix_sha != saved['sha256']:
        raise RuntimeError('PRESERVED_PREFIX_IDENTITY_CHANGED')
    suffix.parent.mkdir(exist_ok=False)
    headers, stats, errors = [suffix.parent / n for n in ['HTTP_headers_private.txt', 'curl_status_private.json', 'curl_stderr_private.txt']]
    with suffix.open('xb'):
        pass
    command = ['/usr/bin/curl', '--disable', '--fail', '--silent', '--show-error',
               '--suppress-connect-headers', '--proto', '=https', '--connect-timeout', '30',
               '--max-time', str(int(deadline-time.monotonic())),
               '--max-filesize', str(limits['maximum_transfer_bytes']),
               '--header', 'Accept-Encoding: identity', '--header', 'If-Match: "'+source['etag']+'"',
               '--range', str(saved['bytes'])+'-'+str(source['bytes']-1),
               '--dump-header', str(headers), '--output', str(suffix), '--write-out', '%{json}', source['url']]
    stop, peak, response, full_md5, full_sha = None, 0, None, None, None
    def guard():
        if time.monotonic() > deadline:
            raise RuntimeError('SHARED_RESUME_DEADLINE')
        if shutil.disk_usage(ROOT).free < floor:
            raise RuntimeError('RESUME_INTERNAL_FLOOR')
        if sum(p.stat().st_size for p in [prefix, suffix, assembled] if p.exists()) > limits['maximum_owned_staging_bytes']:
            raise RuntimeError('RESUME_STAGING_SIZE_LIMIT')
    with stats.open('xb') as sout, errors.open('xb') as serr:
        process = subprocess.Popen(command, stdout=sout, stderr=serr)
        try:
            while process.poll() is None:
                guard()
                rss = subprocess.run(['ps', '-p', str(process.pid), '-o', 'rss='], capture_output=True, text=True).stdout.strip()
                if rss:
                    peak = max(peak, int(rss)*1024)
                if peak > limits['curl_max_rss_bytes']:
                    raise RuntimeError('RESUME_CURL_RSS_LIMIT')
                if headers.exists():
                    response = range_header_check(headers.read_text(), source, saved['bytes'])
                    if response is not None and not response['valid']:
                        raise RuntimeError('WRONG_RANGE_RESPONSE_PRESERVED_NO_ASSEMBLY')
                time.sleep(2)
        except Exception as error:
            stop = str(error)
        finally:
            if process.poll() is None:
                process.terminate()
                try:
                    process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    process.kill(); process.wait()
    try:
        transport = json.loads(stats.read_text())
    except (ValueError, OSError):
        transport = {}
    response = range_header_check(headers.read_text() if headers.exists() else '', source, saved['bytes'])
    verified = False
    try:
        guard()
        if (stop or process.returncode != 0 or transport.get('http_code') != 206
                or not response or not response['valid']
                or suffix.stat().st_size != source['bytes']-saved['bytes']):
            raise RuntimeError(stop or 'RANGE_TRANSFER_NOT_COMPLETE')
        md5, sha = hashlib.md5(), hashlib.sha256()
        prefix_after = hashlib.sha256()
        with assembled.open('xb') as dest:
            for part in [prefix, suffix]:
                with part.open('rb') as source_file:
                    for block in iter(lambda: source_file.read(1024**2), b''):
                        guard(); dest.write(block); md5.update(block); sha.update(block)
                        if part == prefix:
                            prefix_after.update(block)
            dest.flush(); os.fsync(dest.fileno())
        full_md5, full_sha = md5.hexdigest(), sha.hexdigest()
        if (prefix.stat().st_size != saved['bytes'] or prefix_after.hexdigest() != saved['sha256']
                or assembled.stat().st_size != source['bytes'] or full_md5 != source['md5_from_yaml']
                or plan_path.read_bytes() != plan_bytes):
            raise RuntimeError('ASSEMBLED_SOURCE_OR_PREFIX_IDENTITY_FAILED')
        guard()
        publish = subprocess.run(['mv', '-n', str(assembled), str(out)], capture_output=True, text=True)
        if publish.returncode != 0 or assembled.exists() or not out.exists():
            raise RuntimeError('NON_CLOBBERING_PUBLICATION_NOT_COMPLETED')
        verified = True
    except Exception as error:
        stop = str(error)
    receipt = {'completed_utc': datetime.now(timezone.utc).isoformat(),
               'status': 'EXACT_SOURCE_ACQUISITION_PASS' if verified else 'FAILED_RESUME_CHECKPOINTS_PRESERVED',
               'plan_sha256': plan_sha, 'script_sha256': script_sha, 'url': source['url'],
               'path': str(out if verified else suffix), 'prefix_path': str(prefix),
               'prefix_bytes': saved['bytes'], 'prefix_sha256_at_start': prefix_sha,
               'suffix_bytes': suffix.stat().st_size, 'expected_bytes': source['bytes'],
               'actual_bytes': out.stat().st_size if verified else None,
               'expected_md5': source['md5_from_yaml'], 'actual_md5': full_md5, 'actual_sha256': full_sha,
               'exit_code': process.returncode, 'http_code': transport.get('http_code'),
               'range_response': response, 'stop_reason': stop,
               'elapsed_seconds': time.monotonic()-started, 'curl_peak_observed_rss_bytes': peak,
               'internal_free_after_bytes': shutil.disk_usage(ROOT).free,
               'native_3GiB_guard_unchanged': True, 'source_only_acquisition': True,
               'original_sources_modified': False, 'h2_or_rg_runs': 0,
               'preprocessing_admitted': False, 'downstream_source_QC_status': 'NOT_TESTED'}
    persist_receipt(receipt, receipt_path, ssd_receipt_path)
    print(json.dumps(receipt), flush=True)
    if not verified:
        raise RuntimeError('RESUME_NOT_VERIFIED_CHECKPOINTS_PRESERVED')


if __name__ == '__main__':
    if len(sys.argv) == 3 and sys.argv[1] == '--resume-plan':
        resume_main(sys.argv[2])
    elif len(sys.argv) == 1:
        main()
    else:
        raise SystemExit('Usage: acquire_gera_raw.py [--resume-plan PATH]')
