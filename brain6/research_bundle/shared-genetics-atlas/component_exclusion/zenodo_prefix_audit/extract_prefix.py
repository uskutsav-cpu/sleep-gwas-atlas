import zlib,csv,io,json,hashlib,pathlib
p=pathlib.Path(__file__).resolve().parent
raw=(p/'exhaustion_first_131072_bytes.bin').read_bytes();receipt=json.loads((p/'prefix_receipt.json').read_text());assert len(raw)==131072 and hashlib.sha256(raw).hexdigest()==receipt['prefix_sha256']
z=zlib.decompressobj(31);s=z.decompress(raw).decode();rows=list(csv.DictReader(io.StringIO(s[:s.rfind('\n')]),delimiter='\t'))
eligible=[r for r in rows if len(r['ALLELE0'])==len(r['ALLELE1'])==1 and set((r['ALLELE0'],r['ALLELE1'])) not in ({'A','T'},{'C','G'}) and r['ID'].startswith('rs')];chosen=[eligible[round(i*(len(eligible)-1)/19)] for i in range(20)]
with (p/'selected_rows.tsv').open() as f: saved=list(csv.DictReader(f,delimiter='\t'))
assert chosen==saved and len(rows)==receipt['complete_rows'] and len(eligible)==receipt['eligible_nonpalindromic_snv_rows'];print('verified',len(rows),'complete rows;',len(chosen),'selected rows match')
