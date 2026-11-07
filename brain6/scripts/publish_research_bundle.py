#!/usr/bin/env python3
"""Copy the local research delivery into a new, portable repository snapshot.

Downloaded source documents, GWAS inputs and large intermediate matrices stay
local. Record their exact identities rather than distributing them implicitly.
"""
import argparse
import csv
import hashlib
from pathlib import Path
import re
import shutil
from urllib.parse import unquote, urlsplit

ap = argparse.ArgumentParser()
ap.add_argument('--outputs', type=Path, required=True)
args = ap.parse_args()
repo = Path(__file__).resolve().parents[2]
source = args.outputs.resolve()/'shared-genetics-atlas'
bundle = repo/'brain6/research_bundle'
target = bundle/'shared-genetics-atlas'
target.mkdir(parents=True, exist_ok=False)
(bundle/'brain6-alphagenome-pilot').symlink_to('../alphagenome_pilot', target_is_directory=True)
excluded = []
retained = []
for f in sorted(source.rglob('*')):
    if not f.is_file():
        continue
    rel = f.relative_to(source)
    reason = None
    if '__pycache__' in rel.parts or f.suffix == '.pyc' or f.name.startswith('._'):
        continue
    if f.suffix in {'.gz','.xlsx','.zip','.pdf','.docx','.npz','.tbi','.bin','.bed'}:
        reason = 'Downloaded source file or regenerable GWAS/LD intermediate; retained in local research delivery'
    if f.name == 'NG00067.v21.csv':
        reason = 'Third-party individual reference data; not redistributed'
    if '.partial.' in f.name:
        reason = 'Transient partial output'
    if reason:
        excluded.append((str(rel),f.stat().st_size,hashlib.sha256(f.read_bytes()).hexdigest(),reason))
    else:
        out = target/rel
        out.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(f,out)
        retained.append(out)
with (target/'external_files.tsv').open('w',newline='') as f:
    w = csv.writer(f,delimiter='\t')
    w.writerow(['original_path_relative_to_atlas','bytes','sha256','reason'])
    w.writerows(excluded)

# Make navigation portable. Keep the original target visible when a linked
# source is intentionally external; do not create dummy evidence files.
outside = []
for f in retained:
    if f.suffix != '.md':
        continue
    def replace(match):
        label, raw = match.group(1),match.group(2)
        raw = raw.strip().strip('<>')
        u = urlsplit(raw)
        if u.scheme or not u.path or raw.startswith('#'):
            return match.group(0)
        dest = Path(unquote(u.path)) if u.path.startswith('/') else f.parent/unquote(u.path)
        try:
            inside = dest.resolve().is_relative_to(repo.resolve())
        except OSError:
            inside = False
        if inside and dest.exists():
            if 'brain6-alphagenome-pilot' in dest.parts:
                # GitHub's file browser does not traverse directory symlinks.
                canonical = Path(__import__('os').path.relpath(dest.resolve(),f.parent)).as_posix()
                if u.fragment:
                    canonical += '#'+u.fragment
                return f'[{label}]({canonical})'
            return match.group(0)
        outside.append((str(f.relative_to(target)),raw))
        guide = Path(__import__('os').path.relpath(target/'PUBLICATION.md',f.parent)).as_posix()
        return f'[{label}]({guide}#external-dependencies) (external target: `{raw}`)'
    text = f.read_text(errors='replace')
    f.write_text(re.sub(r'(?<!!)\[([^\]\n]*)\]\(([^)\n]+)\)',replace,text))
with (target/'external_links.tsv').open('w',newline='') as f:
    w=csv.writer(f,delimiter='\t')
    w.writerow(['report_path','original_external_or_omitted_target'])
    w.writerows(outside)
print(f'Copied {len(retained)} research files; recorded {len(excluded)} external files and {len(outside)} external links.')
