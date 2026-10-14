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

def digest_prefix_segments(segments, deadline, internal_floor):
    """Hash a contiguous prefix from existing pieces without copying a new prefix."""
    combined_md5, combined_sha = hashlib.md5(), hashlib.sha256()
    identities = []
    for spec in segments:
        path = Path(spec['path'])
        segment_md5, segment_sha = hashlib.md5(), hashlib.sha256()
        size = 0
        with path.open('rb') as stream:
            for block in iter(lambda: stream.read(1024**2), b''):
                if time.monotonic() > deadline:
                    raise RuntimeError('VERIFICATION_WALL_TIME_GUARD')
                if shutil.disk_usage(ROOT).free < internal_floor:
                    raise RuntimeError('VERIFICATION_INTERNAL_FLOOR')
                size += len(block)
                segment_md5.update(block); segment_sha.update(block)
                combined_md5.update(block); combined_sha.update(block)
        if size != spec['bytes']:
            raise RuntimeError('PRESERVED_PREFIX_SEGMENT_LENGTH_CHANGED')
        actual_sha = segment_sha.hexdigest()
        if spec.get('sha256') is not None and actual_sha != spec['sha256']:
            raise RuntimeError('PRESERVED_PREFIX_SEGMENT_IDENTITY_CHANGED')
        identities.append({'path': str(path), 'bytes': size,
                           'md5': segment_md5.hexdigest(), 'sha256': actual_sha})
    return combined_md5.hexdigest(), combined_sha.hexdigest(), identities

def preserved_staging_bytes(files, maximum, reserve):
    """Validate retained files and include macOS AppleDouble companions."""
    seen, total = set(), 0
    for item in files:
        path = Path(item['path'])
        key = str(path)
        if key in seen:
            raise RuntimeError('DUPLICATE_PRESERVED_STAGING_PATH')
        seen.add(key)
        st = path.stat()
        if (not path.is_file() or st.st_size != item['bytes']
                or (item.get('mtime_ns') is not None and st.st_mtime_ns != item['mtime_ns'])):
            raise RuntimeError('PRESERVED_STAGING_INVENTORY_CHANGED')
        total += st.st_size
        sidecar = path.parent / ('._' + path.name)
        if sidecar.exists():
            if not sidecar.is_file():
                raise RuntimeError('INVALID_PRESERVED_STAGING_SIDECAR')
            total += sidecar.stat().st_size
    if total + reserve > maximum:
        raise RuntimeError('PRESERVED_STAGING_EXCEEDS_GLOBAL_CAP')
    return total

def resume_staging_bytes(retained_files, active_files, staging_dir,
                         maximum, metadata_reserve, receipt_name):
    """Count preserved inputs, live data, and every file in the fresh stage."""
    retained_bytes = preserved_staging_bytes(retained_files, maximum, 0)
    retained_metadata_bytes = sum(
        (Path(item['path']).parent / ('._' + Path(item['path']).name)).stat().st_size
        for item in retained_files
        if (Path(item['path']).parent / ('._' + Path(item['path']).name)).exists())
    active_bytes = sum(path.stat().st_size for path in active_files if path.exists())
    actual_metadata_bytes = 0
    if staging_dir.exists():
        allowed_names = {'suffix.part', 'assembled_source.part', 'HTTP_headers_private.txt',
                         'curl_status_private.json', 'curl_stderr_private.txt', receipt_name}
        for path in staging_dir.iterdir():
            name = path.name[2:] if path.name.startswith('._') else path.name
            if name not in allowed_names or not path.is_file():
                raise RuntimeError('UNEXPECTED_RESUME_STAGING_ARTIFACT')
            if path not in active_files or path.name.startswith('._'):
                actual_metadata_bytes += path.stat().st_size
    if retained_metadata_bytes + actual_metadata_bytes > metadata_reserve:
        raise RuntimeError('RESUME_METADATA_RESERVE_EXCEEDED')
    total = retained_bytes + active_bytes + actual_metadata_bytes
    if total > maximum:
        raise RuntimeError('RESUME_STAGING_SIZE_LIMIT')
    return total

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
    prefix_specs = [{'path': spec['path'], 'bytes': spec['bytes'],
                     'sha256': spec.get('sha256'), 'mtime_ns': spec.get('mtime_ns'),
                     'range_start': spec.get('range_start'), 'range_end': spec.get('range_end')}
                    for spec in plan.get('prefix_segments', [saved])]
    prefix_bytes = sum(spec['bytes'] for spec in prefix_specs)
    if 'prefix_segments' in plan:
        position = 0
        for spec in prefix_specs:
            if (spec['range_start'] != position or spec['range_end'] != position + spec['bytes'] - 1):
                raise RuntimeError('PREFIX_SEGMENTS_NOT_EXACTLY_CONTIGUOUS')
            position += spec['bytes']
        if len({spec['path'] for spec in prefix_specs}) != len(prefix_specs):
            raise RuntimeError('DUPLICATE_PREFIX_SEGMENT_PATH')
    plan_sha = hashlib.sha256(plan_bytes).hexdigest()
    script_sha = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    review = json.loads((ROOT / plan['review']).read_text())
    if (review['reviewed_plan_sha256'] != plan_sha or review['reviewed_script_sha256'] != script_sha
            or review['resume_admitted'] is not True):
        raise RuntimeError('EXACT_RESUME_ADMISSION_REQUIRED')
    prefix_paths = [Path(spec['path']) for spec in prefix_specs]
    prefix, suffix = prefix_paths[0], Path(plan['suffix'])
    assembled, out = Path(plan['assembled']), Path(plan['destination'])
    receipt_path = ROOT / plan['receipt']
    ssd_receipt_path = suffix.parent / receipt_path.name
    if any(p.exists() for p in [suffix.parent, assembled, out, receipt_path, ssd_receipt_path]):
        raise RuntimeError('PRIOR_PREFIX_SUFFIX_SOURCE_OR_RECEIPT_PRESERVED')
    started = time.monotonic()
    deadline = started + limits['maximum_seconds']
    floor = limits['internal_emergency_floor_bytes']
    if shutil.disk_usage(ROOT).free < limits['internal_minimum_before_bytes']:
        raise RuntimeError('RESUME_INTERNAL_GUARD')
    if shutil.disk_usage('/Volumes/Extreme SSD').free < limits['ssd_minimum_before_bytes']:
        raise RuntimeError('RESUME_SSD_GUARD')
    retained = limits.get('retained_staging_files', [])
    retained_bytes = preserved_staging_bytes(
        retained, limits['maximum_owned_staging_bytes'], 0)
    if (limits.get('expected_retained_staging_bytes') is not None
            and sum(item['bytes'] for item in retained)
            != limits['expected_retained_staging_bytes']):
        raise RuntimeError('PRESERVED_STAGING_INVENTORY_TOTAL_CHANGED')
    retained_paths = {str(Path(item['path'])) for item in retained}
    if retained and any(str(path) not in retained_paths for path in prefix_paths):
        raise RuntimeError('PREFIX_SEGMENT_MISSING_FROM_STAGING_INVENTORY')
    segment_by_path = {spec['path']: spec for spec in prefix_specs}
    for item in retained:
        if item['path'] in segment_by_path and item['bytes'] != segment_by_path[item['path']]['bytes']:
            raise RuntimeError('PREFIX_STAGING_INVENTORY_LENGTH_MISMATCH')
    prefix_md5, prefix_sha, prefix_identities = digest_prefix_segments(prefix_specs, deadline, floor)
    if prefix_bytes <= 0 or prefix_bytes > source['bytes']:
        raise RuntimeError('PRESERVED_PREFIX_TOTAL_LENGTH_INVALID')
    if len(prefix_specs) == 1 and prefix_sha != saved['sha256']:
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
               '--range', str(prefix_bytes)+'-'+str(source['bytes']-1),
               '--dump-header', str(headers), '--output', str(suffix), '--write-out', '%{json}', source['url']]
    stop, peak, response, full_md5, full_sha = None, 0, None, None, None
    def guard():
        if time.monotonic() > deadline:
            raise RuntimeError('SHARED_RESUME_DEADLINE')
        if shutil.disk_usage(ROOT).free < floor:
            raise RuntimeError('RESUME_INTERNAL_FLOOR')
        active_files = [suffix, assembled, out]
        if not retained:
            active_files.extend(prefix_paths)
        resume_staging_bytes(retained, active_files, suffix.parent,
                             limits['maximum_owned_staging_bytes'],
                             limits.get('reserved_metadata_bytes', 0), receipt_path.name)
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
                    response = range_header_check(headers.read_text(), source, prefix_bytes)
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
    response = range_header_check(headers.read_text() if headers.exists() else '', source, prefix_bytes)
    verified = False
    try:
        guard()
        if (stop or process.returncode != 0 or transport.get('http_code') != 206
                or not response or not response['valid']
                or suffix.stat().st_size != source['bytes']-prefix_bytes):
            raise RuntimeError(stop or 'RANGE_TRANSFER_NOT_COMPLETE')
        md5, sha = hashlib.md5(), hashlib.sha256()
        prefix_after = hashlib.sha256()
        prefix_segment_after = []
        with assembled.open('xb') as dest:
            for index, part in enumerate(prefix_paths + [suffix]):
                part_sha = hashlib.sha256()
                with part.open('rb') as source_file:
                    for block in iter(lambda: source_file.read(1024**2), b''):
                        guard(); dest.write(block); md5.update(block); sha.update(block)
                        if index < len(prefix_paths):
                            prefix_after.update(block)
                            part_sha.update(block)
                if index < len(prefix_paths):
                    prefix_segment_after.append(part_sha.hexdigest())
            dest.flush(); os.fsync(dest.fileno())
        full_md5, full_sha = md5.hexdigest(), sha.hexdigest()
        if (prefix_after.hexdigest() != prefix_sha
                or prefix_segment_after != [item['sha256'] for item in prefix_identities]
                or assembled.stat().st_size != source['bytes'] or full_md5 != source['md5_from_yaml']
                or plan_path.read_bytes() != plan_bytes):
            raise RuntimeError('ASSEMBLED_SOURCE_OR_PREFIX_IDENTITY_FAILED')
        guard()
        publish = subprocess.run(['mv', '-n', str(assembled), str(out)], capture_output=True, text=True)
        if publish.returncode != 0 or assembled.exists() or not out.exists():
            raise RuntimeError('NON_CLOBBERING_PUBLICATION_NOT_COMPLETED')
        guard()
        verified = True
    except Exception as error:
        stop = str(error)
    receipt = {'completed_utc': datetime.now(timezone.utc).isoformat(),
               'status': 'EXACT_SOURCE_ACQUISITION_PASS' if verified else 'FAILED_RESUME_CHECKPOINTS_PRESERVED',
               'plan_sha256': plan_sha, 'script_sha256': script_sha, 'url': source['url'],
               'path': str(out if verified else suffix), 'prefix_path': str(prefix),
               'prefix_bytes': prefix_bytes, 'prefix_sha256_at_start': prefix_sha,
               'prefix_md5_at_start': prefix_md5,
               'prefix_segment_identities': prefix_identities,
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
