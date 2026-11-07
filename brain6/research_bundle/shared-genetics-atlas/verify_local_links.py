#!/usr/bin/env python3
"""Check simple inline Markdown local link targets; anchors and generated HTML aren't validated."""
import pathlib,re,sys,urllib.parse
root=pathlib.Path(__file__).resolve().parent;missing=[];count=0
for path in root.rglob('*.md'):
 text=path.read_text(encoding='utf-8',errors='replace')
 for target in re.findall(r'(?<!!)\[[^\]\n]*\]\(([^)\n]+)\)',text):
  target=target.strip().strip('<>');url=urllib.parse.urlsplit(target)
  if url.scheme or target.startswith('#') or not url.path:continue
  count+=1;name=urllib.parse.unquote(url.path);dest=pathlib.Path(name) if name.startswith('/') else path.parent/name
  if not dest.exists():missing.append((str(path.relative_to(root)),target))
print(f'local inline markdown links: {count} checked, {len(missing)} missing')
for path,target in missing:print(path,target,sep='\t')
sys.exit(bool(missing))
