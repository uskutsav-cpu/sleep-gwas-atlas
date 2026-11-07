import hashlib, io, json, subprocess
from pathlib import Path
import openpyxl
p = Path(__file__).resolve().parent
url = 'https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-025-59383-0/MediaObjects/41467_2025_59383_MOESM3_ESM.xlsx'
data = subprocess.check_output(['curl', '-fsSL', '--max-time', '50', url])
sha = hashlib.sha256(data).hexdigest()
assert sha == '3bb710a7c1871a79936c0ac35ba283865344aca667c4b770718d46d14a8fe8b1'
book = openpyxl.load_workbook(io.BytesIO(data), read_only=True, data_only=True)
sheet = book['Supplementary Data 4']
rows = list(sheet.iter_rows(values_only=True))
headers = list(rows[3][1:12])
selected = [dict(zip(headers, row[1:12])) for row in rows if row[1:3] == ('Insomnia', 'Physical (summary score)')]
assert len(selected) == 3 and {r['Outcome type'] for r in selected} == {'baseline', 'diff', 'log'}
receipt = {'source_url': url, 'sha256': sha, 'sheet': 'Supplementary Data 4', 'rows': selected}
(p/'prior_mr_receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')
print(json.dumps(selected, indent=2))
