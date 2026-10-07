"""Read-only integrity verification: do not rewrite historical evidence."""
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parents[3]
manifest=ROOT/'brain6/translational_psychiatry_research_v1/research_v1/historical_tracked_hashes.json'
files=json.loads(manifest.read_text());fail=[]
for name,expected in files.items():
    path=ROOT/name
    if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest()!=expected:fail.append(name)
print(json.dumps({'inherited_files_checked':len(files),'drift_count':len(fail),'drift':fail}))
if fail:raise SystemExit(1)
