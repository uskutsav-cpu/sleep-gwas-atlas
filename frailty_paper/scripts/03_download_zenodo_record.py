#!/usr/bin/env python3
from __future__ import annotations
import argparse, hashlib, json, os, subprocess
from pathlib import Path
import requests

def md5(path: Path):
    h = hashlib.md5()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(1024*1024), b''): h.update(b)
    return h.hexdigest()

def download(url, dest: Path, expected_size: int, expected_md5: str):
    part = dest.with_name(dest.name + '.part')
    sidecar = part.with_name(part.name + '.json')
    identity = {'url': url, 'expected_size': expected_size, 'expected_md5': expected_md5}

    def quarantine(path: Path):
        quarantine_path = path.with_name(path.name + '.unverified')
        if quarantine_path.exists():
            quarantine_path = path.with_name(path.name + f'.unverified.{path.stat().st_mtime_ns}')
        path.replace(quarantine_path)
        return quarantine_path

    if dest.exists():
        if dest.stat().st_size == expected_size and md5(dest) == expected_md5:
            print('  existing file checksum OK')
            return
        quarantined = quarantine(dest)
        print(f'  moved unverified prior file to {quarantined.name}')
    if part.exists():
        try:
            saved_identity = json.loads(sidecar.read_text(encoding='utf-8'))
        except (OSError, json.JSONDecodeError):
            saved_identity = None
        if saved_identity != identity:
            quarantined = quarantine(part)
            print(f'  moved unrecognized partial to {quarantined.name}')
            sidecar.unlink(missing_ok=True)
    else:
        sidecar.write_text(json.dumps(identity, sort_keys=True) + '\n', encoding='utf-8')
    if not sidecar.exists():
        sidecar.write_text(json.dumps(identity, sort_keys=True) + '\n', encoding='utf-8')
    if part.exists() and part.stat().st_size >= expected_size:
        if part.stat().st_size == expected_size and md5(part) == expected_md5:
            os.replace(part, dest)
            sidecar.unlink(missing_ok=True)
            return
        quarantined = quarantine(part)
        sidecar.unlink(missing_ok=True)
        print(f'  moved unverified partial to {quarantined.name}')
        sidecar.write_text(json.dumps(identity, sort_keys=True) + '\n', encoding='utf-8')
    subprocess.run([
        'curl', '--fail', '--location', '--silent', '--show-error', '--retry', '5', '--retry-delay', '3',
        '--retry-all-errors', '--continue-at', '-', '--header', 'Accept-Encoding: identity', '--output', str(part), url,
    ], check=True)
    actual_size = part.stat().st_size if part.exists() else 0
    if actual_size != expected_size:
        if actual_size > expected_size:
            quarantined = quarantine(part)
            sidecar.unlink(missing_ok=True)
            raise SystemExit(f'Size mismatch for {dest.name}: expected {expected_size}, got {actual_size}; oversized file retained at {quarantined}')
        raise SystemExit(f'Size mismatch for {dest.name}: expected {expected_size}, got {actual_size}; resumable partial retained at {part}')
    actual_md5 = md5(part)
    if actual_md5 != expected_md5:
        quarantined = quarantine(part)
        sidecar.unlink(missing_ok=True)
        raise SystemExit(f'MD5 mismatch for {dest.name}: expected {expected_md5}, got {actual_md5}; file retained at {quarantined}')
    os.replace(part, dest)
    sidecar.unlink(missing_ok=True)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--record', required=True)
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--only', nargs='+', help='Download only these exact filenames from the Zenodo record')
    args=ap.parse_args()
    out=Path(args.outdir); out.mkdir(parents=True, exist_ok=True)
    meta=requests.get(f'https://zenodo.org/api/records/{args.record}', timeout=60)
    meta.raise_for_status(); j=meta.json()
    (out/'zenodo_metadata.json').write_text(json.dumps(j, indent=2), encoding='utf-8')
    items=j.get('files', [])
    if args.only:
        found={item.get('key') or item.get('filename') for item in items}
        unknown=sorted(set(args.only)-found)
        if unknown: raise SystemExit(f'Files not found in Zenodo record {args.record}: {unknown}')
        items=[item for item in items if (item.get('key') or item.get('filename')) in set(args.only)]
    for item in items:
        name=item.get('key') or item.get('filename')
        url=(item.get('links') or {}).get('content') or (item.get('links') or {}).get('self')
        if not name or not url: continue
        dest=out/name
        print('Downloading', name)
        checksum=item.get('checksum','')
        if checksum.startswith('md5:'):
            expected=checksum.split(':',1)[1].lower()
        else:
            raise SystemExit(f'No MD5 registered for {name}; refusing to promote an unverified download')
        download(url, dest, int(item['size']), expected)
        print('  size and md5 OK')
    print('Done:', out)
if __name__=='__main__': main()
