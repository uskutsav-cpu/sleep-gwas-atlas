#!/usr/bin/env python3
"""Reconcile retained PubMed EFetch XML with official ESearch IDs and node types."""
from __future__ import annotations
import argparse,csv,hashlib,json,time,xml.etree.ElementTree as ET
from datetime import datetime,timezone
from pathlib import Path
import requests
BASE='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi'

def xml_records(path):
    ids=[]; counts={'PubmedArticle':0,'PubmedBookArticle':0}
    for _,el in ET.iterparse(path,events=('end',)):
        if el.tag=='PubmedArticle':
            counts['PubmedArticle']+=1; n=el.find('./MedlineCitation/PMID')
        elif el.tag=='PubmedBookArticle':
            counts['PubmedBookArticle']+=1; n=el.find('./BookDocument/PMID')
        else: continue
        if n is not None and n.text: ids.append(n.text.strip())
        el.clear()
    return ids,counts

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--pubmed-dir',default='frailty_paper/review/pubmed'); ap.add_argument('--outdir',default='frailty_paper/review'); a=ap.parse_args()
    root=Path(a.pubmed_dir); out=Path(a.outdir); s=requests.Session(); s.headers['User-Agent']='sleep-frailty-pubmed-reconciliation/1.1'
    segfiles=sorted(root.glob('*/segment_*/query.txt')); rows=[]; queryrows=[]; errors=[]
    for qdir in sorted(p for p in root.iterdir() if p.is_dir() and (p/'query.txt').exists()):
        full=(qdir/'query.txt').read_text(encoding='utf-8').strip()
        try:
            z=s.get(BASE,params={'db':'pubmed','term':full,'retmode':'json','retmax':0},timeout=90); z.raise_for_status(); fullcount=int(z.json()['esearchresult']['count'])
        except Exception as e:
            errors.append({'query_file':str((qdir/'query.txt').relative_to(root)),'error':f'full-query ESearch: {e}'}); fullcount=-1
        queryseg=[p for p in segfiles if p.parent.parent==qdir]
        union_search=set(); union_xml=set(); ssum=xmlsum=articles=books=miss_total=extra_total=0
        for qfile in queryseg:
            rel=str(qfile.relative_to(root)); query=qfile.read_text(encoding='utf-8').strip()
            try:
                z=s.get(BASE,params={'db':'pubmed','term':query,'retmode':'json','retmax':10000},timeout=90); z.raise_for_status(); d=z.json()['esearchresult']; ids=[str(x) for x in d.get('idlist',[])]; count=int(d['count'])
                if count!=len(ids): errors.append({'query_file':rel,'error':f'ESearch count {count} != idlist length {len(ids)}'})
            except Exception as e:
                errors.append({'query_file':rel,'error':f'segment ESearch: {e}'}); continue
            xmls=sorted(qfile.parent.glob('batch_*.xml')); got=[]; c={'PubmedArticle':0,'PubmedBookArticle':0}
            for xp in xmls:
                vals,k=xml_records(xp); got.extend(vals)
                for name in c:c[name]+=k[name]
            exp=set(ids); seen=set(got); missing=sorted(exp-seen); extra=sorted(seen-exp)
            rows.append({'query_file':rel,'esearch_count_now':count,'esearch_ids':len(ids),'xml_files':len(xmls),'xml_article_nodes':c['PubmedArticle'],'xml_book_nodes':c['PubmedBookArticle'],'xml_record_nodes':sum(c.values()),'xml_unique_pmids':len(seen),'missing_pmids':len(missing),'xml_pmids_not_in_esearch':len(extra),'missing_pmid_list':';'.join(missing),'expected_id_sha256':hashlib.sha256('\n'.join(sorted(exp)).encode()).hexdigest()})
            union_search.update(exp); union_xml.update(seen); ssum+=count; xmlsum+=len(got); articles+=c['PubmedArticle']; books+=c['PubmedBookArticle']; miss_total+=len(missing); extra_total+=len(extra)
            time.sleep(.36)
        queryrows.append({'query_id':qdir.name,'full_query_count_now':fullcount,'segment_count_sum':ssum,'segment_union_unique_pmids':len(union_search),'xml_article_nodes':articles,'xml_book_nodes':books,'xml_record_nodes':xmlsum,'xml_union_unique_pmids':len(union_xml),'missing_segment_id_occurrences':miss_total,'extra_segment_id_occurrences':extra_total,'segment_union_minus_full_count':len(union_search)-fullcount if fullcount>=0 else ''})
        time.sleep(.36)
    out.mkdir(parents=True,exist_ok=True)
    def write(name,fields,data):
        with (out/name).open('w',newline='',encoding='utf-8') as f:
            w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(data)
    write('pubmed_retrieval_reconciliation.tsv',['query_file','esearch_count_now','esearch_ids','xml_files','xml_article_nodes','xml_book_nodes','xml_record_nodes','xml_unique_pmids','missing_pmids','xml_pmids_not_in_esearch','missing_pmid_list','expected_id_sha256'],rows)
    write('pubmed_retrieval_reconciliation_by_query.tsv',['query_id','full_query_count_now','segment_count_sum','segment_union_unique_pmids','xml_article_nodes','xml_book_nodes','xml_record_nodes','xml_union_unique_pmids','missing_segment_id_occurrences','extra_segment_id_occurrences','segment_union_minus_full_count'],queryrows)
    write('pubmed_retrieval_reconciliation_errors.tsv',['query_file','error'],errors)
    summary={'checked_at_utc':datetime.now(timezone.utc).isoformat(),'source':'NCBI E-utilities ESearch','endpoint':BASE,'segments_checked':len(rows),'query_summaries':len(queryrows),'query_errors':len(errors),'missing_segment_id_occurrences':sum(int(x['missing_pmids']) for x in rows),'extra_segment_id_occurrences':sum(int(x['xml_pmids_not_in_esearch']) for x in rows),'pubmed_article_nodes':sum(x['xml_article_nodes'] for x in queryrows),'pubmed_book_article_nodes':sum(x['xml_book_nodes'] for x in queryrows),'command':'python frailty_paper/scripts/09_reconcile_pubmed_retrieval.py'}
    (out/'pubmed_retrieval_reconciliation.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
