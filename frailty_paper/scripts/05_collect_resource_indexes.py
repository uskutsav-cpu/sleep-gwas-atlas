#!/usr/bin/env python3
from __future__ import annotations
import argparse, json
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import urljoin
import subprocess

EQTL_META = 'https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/data_tables/dataset_metadata_r7.tsv'
BRAINSCOPE_KEY = 'https://brainscope.gersteinlab.org/key_resource_files.html'
BRAINSCOPE_EQTL = 'https://brainscope.gersteinlab.org/output-sig-eQTL.html'

def save_url(url, path):
    result=subprocess.run(['curl','--fail','--location','--silent','--show-error','--max-time','90',url],check=True,capture_output=True)
    path.write_bytes(result.stdout)

def collect_links(page, out_json):
    result=subprocess.run(['curl','--fail','--location','--silent','--show-error','--max-time','90',page],check=True,capture_output=True,text=True)
    html=result.stdout
    rows=[]
    class Links(HTMLParser):
        def __init__(self):
            super().__init__(); self.rows=[]; self.href=None; self.label=[]
        def handle_starttag(self, tag, attrs):
            if tag=='a':
                self.href=dict(attrs).get('href'); self.label=[]
        def handle_data(self, data):
            if self.href is not None: self.label.append(data.strip())
        def handle_endtag(self, tag):
            if tag=='a' and self.href is not None:
                href=urljoin(page,self.href); label=' '.join(x for x in self.label if x)
                if href.startswith('http'): self.rows.append({'label':label,'url':href})
                self.href=None; self.label=[]
    parser=Links(); parser.feed(html); rows=parser.rows
    out_json.write_text(json.dumps(rows, indent=2), encoding='utf-8')
    return rows

def maybe_download_by_name(rows, names, outdir):
    outdir.mkdir(parents=True, exist_ok=True)
    for name in names:
        hits=[x for x in rows if name.lower() in (x['label']+' '+x['url']).lower()]
        if not hits:
            print('No exact link found for', name); continue
        url=hits[0]['url']; dest=outdir/name
        if dest.exists() and dest.stat().st_size>0:
            print('Exists:', dest); continue
        print('Downloading', name, 'from', url)
        subprocess.run(['curl','--fail','--location','--silent','--show-error','--max-time','120','--output',str(dest),url],check=True)

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--root', required=True); ap.add_argument('--download-brainscope-key', action='store_true'); ap.add_argument('--download-brainscope-sample-metadata', action='store_true'); args=ap.parse_args()
    root=Path(args.root)
    qtl=root/'data/qtl'; sc=root/'data/single_cell'; qtl.mkdir(parents=True, exist_ok=True); sc.mkdir(parents=True, exist_ok=True)
    if args.download_brainscope_sample_metadata:
        link_index=sc/'brainscope_key_resource_links.json'
        if not link_index.is_file(): raise FileNotFoundError(f'Missing cached BrainSCOPE link index: {link_index}')
        key_links=json.loads(link_index.read_text(encoding='utf-8'))
        maybe_download_by_name(key_links, ['PEC2_sample_metadata.txt'], sc/'brainscope_key_files')
        print('Downloaded only the public BrainSCOPE sample metadata; no QTL archive requested.')
        return
    failures=[]
    try:
        save_url(EQTL_META, qtl/'eqtl_catalogue_dataset_metadata_r7.tsv')
        print('Saved eQTL Catalogue metadata index')
    except (subprocess.SubprocessError, OSError) as exc:
        failures.append(f'eQTL Catalogue metadata: {exc}')
        print(f'FAILED eQTL Catalogue metadata: {exc}')
    try:
        key_links=collect_links(BRAINSCOPE_KEY, sc/'brainscope_key_resource_links.json')
        print(f'Saved {len(key_links)} BrainSCOPE key-resource links')
    except (subprocess.SubprocessError, OSError) as exc:
        failures.append(f'BrainSCOPE key-resource index: {exc}')
        key_links=[]
        print(f'FAILED BrainSCOPE key-resource index: {exc}')
    try:
        eqtl_links=collect_links(BRAINSCOPE_EQTL, sc/'brainscope_eqtl_links.json')
        print(f'Saved {len(eqtl_links)} BrainSCOPE eQTL links')
    except (subprocess.SubprocessError, OSError) as exc:
        failures.append(f'BrainSCOPE eQTL index: {exc}')
        eqtl_links=[]
        print(f'FAILED BrainSCOPE eQTL index: {exc}')
    if args.download_brainscope_key:
        maybe_download_by_name(key_links, ['PEC2_sample_metadata.txt'], sc/'brainscope_key_files')
        maybe_download_by_name(eqtl_links, ['sig_QTLs.zip'], sc/'brainscope_key_files')
    print('Resource indexes collected. Full QTL matrices are intentionally NOT bulk-downloaded; query/download them after loci are frozen.')
    if failures:
        raise SystemExit(f'{len(failures)} source index request(s) failed; successful indexes were retained.')
if __name__=='__main__': main()
