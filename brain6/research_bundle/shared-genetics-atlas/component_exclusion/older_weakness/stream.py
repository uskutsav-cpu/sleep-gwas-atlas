import csv, gzip, hashlib, json, math, os, subprocess, zlib
from pathlib import Path

out = Path(__file__).resolve().parent
inputs = Path('/Volumes/Brain6Work/physical-function-source-sensitivity-20261028')
common = {}
with gzip.open(inputs / 'questionnaire.sumstats.gz', 'rt') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        common[row['SNP']] = frozenset((row['A1'], row['A2']))
assert len(common) == 1020637
ref = {}
with gzip.open('/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/reference/hm3_grch37_variant_map.tsv.gz', 'rt') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        if row['SNP'] in common:
            assert row['SNP'] not in ref
            ref[row['SNP']] = (row['CHR'], row['BP'], frozenset((row['A1'], row['A2'])))
url = 'https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90007001-GCST90008000/GCST90007526/GCST90007526_buildGRCh37.tsv.gz'
proc = subprocess.Popen(['curl', '-fsSL', '--retry', '2', url], stdout=subprocess.PIPE)
h = hashlib.md5(); sha = hashlib.sha256(); dec = zlib.decompressobj(31); buf = b''
rows = kept = bad = duplicates = size = 0; head = None; seen = set(); done = False
neff = 4 * 48596 * 207927 / 256523
print('SNP\tA1\tA2\tN\tZ', flush=True)
def process(line):
    global rows, kept, bad, duplicates, head
    cols = line.decode().rstrip('\r').split('\t')
    if head is None:
        head = cols
        assert len(head) == 17 and head[:6] == ['variant_id','chromosome','base_pair_location','effect_allele','other_allele','effect_allele_frequency']
        return
    rows += 1
    if len(cols) != len(head): raise ValueError('invalid row width')
    d = dict(zip(head, cols)); rs = d['variant_id']; a = d['effect_allele']; b = d['other_allele']; pair = frozenset((a,b)); v = ref.get(rs)
    if v is None or (d['chromosome'],d['base_pair_location'],pair) != v or pair != common[rs] or pair in (frozenset('AT'),frozenset('CG')): return
    try:
        beta, se, p, freq = map(float, (d['beta'],d['standard_error'],d['p_value'],d['effect_allele_frequency']))
        assert all(map(math.isfinite,(beta,se,p,freq))) and se > 0 and 0 <= p <= 1 and 0 < freq < 1
    except (ValueError, AssertionError, OverflowError): bad += 1; return
    if rs in seen: duplicates += 1; return
    seen.add(rs); kept += 1
    print(f'{rs}\t{a}\t{b}\t{neff:.6f}\t{beta/se:.12g}')
while chunk := proc.stdout.read(2**20):
    h.update(chunk); sha.update(chunk); size += len(chunk)
    data = dec.decompress(chunk); done = dec.eof
    if dec.unused_data: raise ValueError('unexpected additional gzip member')
    buf += data; lines = buf.split(b'\n'); buf = lines.pop()
    for line in lines:
        if line: process(line)
if buf: process(buf)
assert proc.wait() == 0 and done and h.hexdigest() == '10ed84976ef644b1ac7ceb06cd34e2dc' and duplicates == 0
(out / ('receipt_' + os.environ['STREAM_TAG'] + '.json')).write_text(json.dumps({'url':url,'bytes':size,'md5':h.hexdigest(),'sha256':sha.hexdigest(),'rows':rows,'retained':kept,'invalid_matched':bad,'duplicate_eligible':duplicates,'working_effective_n':neff},indent=2)+'\n')
