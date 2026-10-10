"""Reference-only BED decoding/QC and bounded Gram controls; no estimator."""
import hashlib
import math
import os
from pathlib import Path
import re
import struct
import zlib

NSAMPLE = 503
BYTES_PER_SNP = 126
CHUNK_MAX = 2048
BED_HEADER = b'\x6c\x1b\x01'


def regular(path):
    path = Path(path)
    if not path.is_file() or path.is_symlink() or any(p.is_symlink() for p in path.parents):
        raise RuntimeError('SOURCE_REGULAR_NO_SYMLINK_REQUIRED: ' + str(path))
    return path


def hashes(path, guard=lambda: None):
    md5, sha = hashlib.md5(), hashlib.sha256()
    with regular(path).open('rb') as f:
        for b in iter(lambda: f.read(65536), b''):
            guard(); md5.update(b); sha.update(b)
    return {'bytes': Path(path).stat().st_size, 'md5': md5.hexdigest(), 'sha256': sha.hexdigest()}


def decode_a1(payload, np, n=NSAMPLE):
    if n != NSAMPLE or len(payload) % BYTES_PER_SNP:
        raise RuntimeError('EXACT503_COMPLETE_BED_ROWS_REQUIRED')
    m = len(payload)//BYTES_PER_SNP
    if not 0 < m <= CHUNK_MAX:
        raise RuntimeError('DECODE_CHUNK_OUTSIDE1_TO2048')
    raw = np.frombuffer(payload, dtype=np.uint8).reshape(m, BYTES_PER_SNP)
    if np.any(raw[:, -1] & 0xC0):
        raise RuntimeError('NONZERO503_PADDING_BITS')
    code = ((raw[:, :, None] >> np.array([0, 2, 4, 6], dtype=np.uint8)) & 3).reshape(m, 504)[:, :503]
    return np.array([2, -1, 1, 0], dtype=np.int8)[code].T


def qc_column(dose, np):
    if dose.shape != (503,) or not np.all(np.isin(dose, [-1, 0, 1, 2])):
        raise RuntimeError('INVALID503_DISCRETE_DOSAGE')
    observed = dose[dose >= 0].astype(np.int64)
    n = int(len(observed)); miss = 503-n
    total = int(observed.sum()); total2 = int((observed*observed).sum())
    numerator = n*total2-total*total
    if n == 0: reason = 'ALL_MISSING'
    elif n < 2: reason = 'FEWER_THAN_TWO_OBSERVED'
    elif numerator <= 0: reason = 'EXACT_ZERO_VARIANCE'
    elif miss > 25: reason = 'MISSING_GT25_OF503'
    else: reason = None
    f = total/(2*n) if n else None
    return {'observed_count': n, 'missing_count': miss, 'sum_dose': total, 'sum_dose2': total2,
            'integer_variance_numerator': numerator, 'A1_frequency': f,
            'MAF': min(f, 1-f) if f is not None else None, 'exclusion': reason}


def standardized(dose, sign, np):
    if sign not in (-1, 1): raise RuntimeError('EXACT_ORDERED_OR_SWAPPED_SIGN_REQUIRED')
    q = qc_column(dose, np)
    if q['exclusion'] is not None: raise RuntimeError('INELIGIBLE_REFERENCE_COLUMN: '+q['exclusion'])
    mean = q['sum_dose']/q['observed_count']
    x = dose.astype(np.float64); x[x < 0] = mean; x -= mean
    scale = np.sqrt(np.dot(x, x)/503)
    if not np.isfinite(scale) or scale <= 0: raise RuntimeError('INVALID_UNREPAIRED_REFERENCE_SCALE')
    x = sign*x/scale
    if not np.all(np.isfinite(x)) or abs(float(x.mean())) > 1e-10 or abs(float(np.dot(x, x)/503)-1) > 1e-10:
        raise RuntimeError('REFERENCE_COLUMN_MEAN_OR_NORM_FAILED_NO_REPAIR')
    return x


def allele_sign(target, source):
    a, b = target
    if a not in 'ACGT' or b not in 'ACGT' or len(a) != 1 or len(b) != 1 or a == b:
        return None, 'INVALID_TARGET_ALLELES'
    if {a, b} in ({'A', 'T'}, {'C', 'G'}): return None, 'PALINDROMIC'
    if tuple(source) == tuple(target): return 1, None
    if tuple(source) == tuple(reversed(target)): return -1, None
    return None, 'ALLELE_CONFLICT_NO_COMPLEMENT_RESCUE'


def frozen_control_indices(n):
    if n < 192: raise RuntimeError('G5_INCOMPLETE_FEWER_THAN192_NO_WINDOW_REPLACEMENT')
    first = list(range(64)); middle = list(range((n-64)//2, (n-64)//2+64)); last = list(range(n-64, n))
    if len(set(first+middle+last)) != 192: raise RuntimeError('G5_OVERLAPPING_PRIMARY_WINDOWS')
    scatter = []
    for begin, end in [(0, n//3), (n//3, 2*n//3), (2*n//3, n)]:
        if end-begin < 32: raise RuntimeError('G5_INCOMPLETE_SCATTER')
        scatter.append([begin+(i*(end-begin-1))//31 for i in range(32)])
    if len(set(sum(scatter, []))) != 96: raise RuntimeError('G5_SCATTER_DUPLICATES')
    return {'first64': first, 'middle64': middle, 'last64': last,
            'early32_scatter': scatter[0], 'middle32_scatter': scatter[1], 'late32_scatter': scatter[2]}


def operator(rows, read_g, v, np, chunk=2048, reverse=False, guard=lambda: None):
    if chunk < 1 or chunk > CHUNK_MAX or len(rows) != len(v): raise RuntimeError('OPERATOR_SHAPE_OR_CHUNK_FAILED')
    groups = [list(range(i, min(i+chunk, len(rows)))) for i in range(0, len(rows), chunk)]
    if reverse: groups.reverse()
    gv = np.zeros(503, dtype=np.float64)
    for ix in groups:
        guard(); g = read_g([rows[i] for i in ix]); gv += g @ v[ix]
    result = np.empty(len(rows), dtype=np.float64)
    for ix in groups:
        guard(); g = read_g([rows[i] for i in ix]); result[ix] = g.T @ (gv/503)
    if not np.all(np.isfinite(result)): raise RuntimeError('NONFINITE_OPERATOR_NO_REPAIR')
    return result, gv


def close_array(actual, expected, np):
    return actual.shape == expected.shape and bool(np.all(np.isfinite(actual))) and bool(np.all(np.isfinite(expected))) and bool(np.all(np.abs(actual-expected) <= 1e-10+1e-8*np.abs(expected)))


class StrictGzip:
    """Exactly one verified gzip member; bounded buffering; reject trailing bytes."""
    def __init__(self, path, cap, guard):
        self.file = regular(path).open('rb'); self.z = zlib.decompressobj(16+zlib.MAX_WBITS)
        self.buffer = b''; self.pending = b''; self.done = False; self.total = 0; self.cap = cap; self.guard = guard

    def read(self, n):
        if not 0 <= n <= 65536: raise RuntimeError('BOUNDED_GZIP_READ_REQUIRED')
        while len(self.buffer) < n and not self.done:
            self.guard()
            data = self.pending or self.file.read(65536); self.pending = b''
            if not data:
                if not self.z.eof: raise RuntimeError('GZIP_TRUNCATION_OR_CRC_EOF')
                self.done = True; break
            result = self.z.decompress(data, max(1, n-len(self.buffer)))
            self.pending = self.z.unconsumed_tail
            self.total += len(result)
            if self.total > self.cap: raise RuntimeError('TAR_DECOMPRESSED_CAP')
            self.buffer += result
            if self.z.eof:
                if self.z.unused_data or self.pending or self.file.read(1): raise RuntimeError('GZIP_EXTRA_MEMBER_OR_TRAILING_PAYLOAD')
                self.done = True
        result, self.buffer = self.buffer[:n], self.buffer[n:]
        return result

    def close(self): self.file.close()


def _octal(field):
    text = field.rstrip(b'\0 ').lstrip(b' ')
    if text and not re.fullmatch(b'[0-7]+', text): raise RuntimeError('TAR_NON_OCTAL_OR_BASE256_METADATA')
    return int(text or b'0', 8)


def extract_strict66(archive, dest, cap, guard=lambda: None):
    """Manual owned extraction; no tarfile extended metadata or extractall."""
    dest = Path(dest)
    if dest.exists() or dest.is_symlink() or any(p.is_symlink() for p in dest.parents): raise RuntimeError('FRESH_OWNED_EXTRACTION_REQUIRED')
    dest.mkdir(); reader = StrictGzip(archive, cap+1024*1024, guard)
    expected = {f'1000G.EUR.QC.{c}.{s}' for c in range(1, 23) for s in ['bed', 'bim', 'fam']}
    seen = set(); prefixes = set(); files = {}; total = 0; zero = 0; directory_names = set()
    try:
        while True:
            h = reader.read(512)
            if len(h) != 512: raise RuntimeError('TAR_SHORT_HEADER_OR_MISSING_TERMINATOR')
            if h == b'\0'*512:
                zero += 1
                if zero == 2: break
                continue
            if zero: raise RuntimeError('TAR_SINGLE_ZERO_THEN_PAYLOAD')
            if sum(h[:148])+8*32+sum(h[156:]) != _octal(h[148:156]): raise RuntimeError('TAR_HEADER_CHECKSUM')
            name = h[:100].split(b'\0')[0].decode('ascii'); prefix = h[345:500].split(b'\0')[0].decode('ascii')
            if prefix: name = prefix+'/'+name
            if name in seen: raise RuntimeError('TAR_DUPLICATE_NAME')
            seen.add(name)
            if name.startswith('/') or '\\' in name or '..' in name.split('/') or '\0' in name: raise RuntimeError('TAR_PATH_REJECTED')
            size = _octal(h[124:136]); kind = h[156:157]
            if kind == b'5':
                if size or name not in ['.', './', '1000G_EUR_Phase3_plink', '1000G_EUR_Phase3_plink/']: raise RuntimeError('TAR_UNEXPECTED_DIRECTORY')
                directory_names.add(name); continue
            if kind not in [b'0', b'\0'] or h[157:257].strip(b'\0'): raise RuntimeError('TAR_LINK_SPECIAL_SPARSE_OR_EXTENDED_METADATA')
            family = '1000G_EUR_Phase3_plink/' if name.startswith('1000G_EUR_Phase3_plink/') else ''
            base = name[len(family):]
            if base not in expected or '/' in base or base in files: raise RuntimeError('TAR_UNEXPECTED_OR_DUPLICATE_FILE')
            prefixes.add(family)
            if len(prefixes) != 1 or len(files) >= 66 or size <= 0 or total+size > cap: raise RuntimeError('TAR_PREFIX_OR_EXTRACTION_CAP')
            target = dest/base; digest = hashlib.sha256(); left = size
            with target.open('xb') as out:
                while left:
                    guard(); b = reader.read(min(65536, left))
                    if not b: raise RuntimeError('TAR_MEMBER_TRUNCATED')
                    out.write(b); digest.update(b); left -= len(b)
                out.flush(); os.fsync(out.fileno())
            padding = (-size) % 512
            if padding and reader.read(padding) != b'\0'*padding: raise RuntimeError('TAR_MEMBER_NONZERO_PADDING')
            total += size; files[base] = {'path': str(target), 'bytes': size, 'sha256': digest.hexdigest(), 'archive_member': name}
        while True:
            b = reader.read(65536)
            if not b: break
            if any(b): raise RuntimeError('TAR_NONZERO_TRAILING_PAYLOAD')
        if set(files) != expected: raise RuntimeError('EXACT66_SOURCE_MEMBERS_REQUIRED')
        if any(n.startswith('1000G_EUR_Phase3_plink') for n in directory_names) and prefixes != {'1000G_EUR_Phase3_plink/'}: raise RuntimeError('TAR_DIRECTORY_PREFIX_DIFFERS')
        fd = os.open(dest, os.O_RDONLY)
        try: os.fsync(fd)
        finally: os.close(fd)
        return {'files': files, 'prefix': next(iter(prefixes)), 'bytes': total, 'gzip_single_member_verified': True}
    finally: reader.close()
