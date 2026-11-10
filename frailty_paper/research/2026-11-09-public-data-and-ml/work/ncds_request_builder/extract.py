import csv, json
from collections import defaultdict
from pathlib import Path
import openpyxl

root = Path(__file__).resolve().parents[2]
inputs = defaultdict(list)
with (root / 'outputs/ncds_age50_frailty_candidate_questionnaire_map.csv').open() as f:
    for row in csv.DictReader(f):
        inputs[row['candidate_variable'].lower()].append((row['published_deficit'], row['match_status']))
groups = {
    'Type 2 diabetes: type not resolved': ['n8khpb03'],
    'Stomach problems: screening and subtypes': ['n8khpb14', 'n8sbgb01', 'n8sbgb02', 'n8sbgb06', 'n8sbgb08', 'n8sbgb09', 'n8sbgb11'],
    'Hearing problems: screening and subtypes': ['n8khpb07', 'n8hprb01', 'n8hprb02', 'n8hprb04', 'n8hprb05', 'n8hprb07', 'n8hprb09'],
    'Eyesight problems: screening and subtypes': ['n8khpb08', 'n8eyep07', 'n8eyep08', 'n8eyep09', 'n8eyep10', 'n8eyep13'],
    'Depression: medical contact and persistence': ['n8mhprb1', 'n8mhstl1'],
    'Anxiety: medical contact and persistence': ['n8mhprb2', 'n8mhstl2'],
    'Sleep maintenance outcome candidate': ['n8scq14'],
}
for label, names in groups.items():
    for name in names:
        inputs[name].append((label, 'candidate; routing and original derivation unverified'))
wb = openpyxl.load_workbook(root / 'work/ncds_cls_dictionary_current.xlsx', read_only=True, data_only=True)
ws = wb['NCDS_Safeguarded_DataDict']
found = defaultdict(list)
for row in ws.iter_rows(min_row=2, values_only=True):
    if row[10] is not None and str(row[10]).lower() in inputs:
        found[str(row[10]).lower()].append(row)
records = []
for var, mappings in inputs.items():
    rows = found.get(var, [])
    if not rows:
        raise ValueError(f'{var}: absent from dictionary')
    for r in rows:
        records.append([str(r[10]), str(r[3] or ''), str(r[4] or ''), str(r[8] or ''), '; '.join(x[0] for x in mappings), '; '.join(x[1] for x in mappings) + ('; duplicated name: select correct dataset with CLS' if len(rows) > 1 else '')])
out = {'source': 'https://cls.ucl.ac.uk/wp-content/uploads/2017/02/NCDS_Longitudinal_Data_Dictionary.xlsx', 'fields': records}
(root / 'work/ncds_request_builder/fields.json').write_text(json.dumps(out, indent=2) + '\n')
print(f'{len(inputs)} unique candidate names; {len(records)} matching dictionary rows')
