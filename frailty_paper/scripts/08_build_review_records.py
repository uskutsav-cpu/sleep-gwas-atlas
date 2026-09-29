#!/usr/bin/env python3
"""Build auditable PRISMA screening files from retained PubMed XML records."""
from __future__ import annotations
import argparse, csv, hashlib, json, re, unicodedata, xml.etree.ElementTree as ET
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

FIELDS = ['record_id','record_type','pmid','doi','normalized_title','title','abstract','year','journal','query_id','source_file','source_database','source_record_id','authors','import_status','duplicate_group','retained_record_id','dedup_rule']
EXTRACTION = ['author','year','PMID','DOI','country','cohort','study_design','sample_size','age','sex','ancestry','sleep_trait','sleep_measure','frailty_definition','frailty_instrument','direction','followup','effect_type','effect','lower_ci','upper_ci','p_value','adjustments','genetic_method','locus','gene','tissue','cell_type','main_conclusion','risk_of_bias','notes']

def sha(path: Path) -> str:
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
    return h.hexdigest()

def pubmed_xml_files(base: Path) -> list[Path]:
    """Return XML records while ignoring macOS AppleDouble sidecars."""
    return sorted(
        path for path in base.rglob('*.xml')
        if path.is_file() and not path.name.startswith('._')
    )

def pubmed_source_path(path: Path, base: Path, source_prefix: Path) -> str:
    """Keep source locators stable when the cache directory is a symlink."""
    return (source_prefix / path.relative_to(base)).as_posix()

def norm(s: str) -> str:
    s=unicodedata.normalize('NFKD',s).encode('ascii','ignore').decode().lower()
    return ' '.join(re.findall(r'[a-z0-9]+',s))

def txt(node) -> str:
    return ' '.join(''.join(node.itertext()).split()) if node is not None else ''

def get_year(article, record_type: str) -> str:
    if record_type=='PubmedArticle':
        for node in article.findall('.//PubDate/Year'):
            if node.text: return node.text.strip()
        med=article.find('.//PubDate/MedlineDate')
    else:
        for node in article.findall('./BookDocument/Book/PubDate/Year'):
            if node.text: return node.text.strip()
        med=article.find('./BookDocument/Book/PubDate/MedlineDate')
    m=re.search(r'\b(19|20)\d{2}\b',txt(med)) if med is not None else None
    return m.group(0) if m else ''

def parse_article(article, path: Path, base: Path, idx: int, source_prefix: Path | None = None):
    record_type=article.tag
    if record_type=='PubmedArticle':
        doc=article.find('./MedlineCitation')
        content=doc.find('./Article') if doc is not None else None
        pmid=txt(doc.find('./PMID')) if doc is not None else ''
        journal=txt(content.find('./Journal/Title')) if content is not None else ''
        abstract_nodes=content.findall('./Abstract/AbstractText') if content is not None else []
        title=txt(content.find('./ArticleTitle')) if content is not None else ''
    else:
        doc=article.find('./BookDocument')
        content=doc.find('./Book') if doc is not None else None
        pmid=txt(doc.find('./PMID')) if doc is not None else ''
        journal=(txt(content.find('./CollectionTitle')) or txt(content.find('./BookTitle'))) if content is not None else ''
        abstract_nodes=doc.findall('./Abstract/AbstractText') if doc is not None else []
        title=(txt(doc.find('./ArticleTitle')) or txt(content.find('./BookTitle'))) if doc is not None else ''
    doi=''
    for node in article.findall('.//ArticleIdList/ArticleId'):
        if node.attrib.get('IdType','').lower()=='doi': doi=txt(node); break
    abst=' '.join(txt(x) for x in abstract_nodes if txt(x))
    parts=path.relative_to(base).parts
    query=parts[0]
    prefix=source_prefix if source_prefix is not None else Path(base.name)
    rel=pubmed_source_path(path,base,prefix)
    occurrence=hashlib.sha256(f'{rel}\n{idx}'.encode()).hexdigest()[:12]
    authors='; '.join(txt(x) for x in article.findall('.//AuthorList/Author') if txt(x))
    record_id=f'{query}:{pmid or "NO-PMID"}:{occurrence}'
    return {'record_id':record_id,'record_type':'journal_article' if record_type=='PubmedArticle' else 'book_article','pmid':pmid,'doi':doi,'normalized_title':norm(title),'title':title,'abstract':abst,'year':get_year(article,record_type),'journal':journal,'query_id':query,'source_file':rel,'source_database':'PubMed','source_record_id':pmid or record_id,'authors':authors,'import_status':'READY' if title else 'MISSING_TITLE'}

def dedup_key(r):
    if r['pmid']: return ('PMID',r['pmid'])
    if r['doi']:
        doi=re.sub(r'^https?://(?:dx\.)?doi\.org/','',r['doi'].strip().lower())
        return ('DOI',doi)
    if r['normalized_title']: return ('TITLE',r['normalized_title'])
    return ('RECORD',r['record_id'])

def dedup_keys(r):
    keys=[]
    if r.get('pmid'): keys.append(('PMID',r['pmid'].strip()))
    if r.get('doi'):
        doi=re.sub(r'^https?://(?:dx\.)?doi\.org/','',r['doi'].strip().lower())
        keys.append(('DOI',doi))
    if r.get('normalized_title'): keys.append(('TITLE',r['normalized_title']))
    if not keys: keys.append(('RECORD',r['record_id']))
    return keys

def group_duplicate_records(records):
    """Deduplicate by PMID, then DOI; use title only without ID conflicts."""
    parent=list(range(len(records)))
    def find(index):
        while parent[index]!=index:
            parent[index]=parent[parent[index]]
            index=parent[index]
        return index
    def union(left,right):
        a,b=find(left),find(right)
        if a!=b: parent[b]=a
    first_for_key={}
    # Group by PMID first.
    for index,record in enumerate(records):
        for key in dedup_keys(record):
            if key[0]!='PMID': continue
            if key in first_for_key: union(index,first_for_key[key])
            else: first_for_key[key]=index
    # DOI can bridge a source lacking PMID, but PubMed metadata sometimes
    # assigns one DOI to records with different titles. Partition such DOI
    # buckets by exact normalized title instead of collapsing them together.
    by_doi=defaultdict(list)
    for index,record in enumerate(records):
        for key in dedup_keys(record):
            if key[0]=='DOI': by_doi[key[1]].append(index)
    for doi,indices in by_doi.items():
        title_groups=defaultdict(list)
        missing_title=[]
        for index in indices:
            title=records[index].get('normalized_title','')
            (title_groups[title] if title else missing_title).append(index)
        for title,group in title_groups.items():
            anchor=group[0]
            for index in group[1:]: union(anchor,index)
        if len(title_groups)==1:
            anchor=next(iter(title_groups.values()))[0]
            for index in missing_title: union(anchor,index)
        elif missing_title:
            # With multiple DOI/title candidates, title-less records remain
            # separate for manual reconciliation rather than guessed matches.
            pass
    # Title fallback is safe only if the normalized-title bucket contains at
    # most one component with any PMID/DOI. Generic titles across distinct
    # papers must remain separate. Identifier-free records can still dedup by
    # exact normalized title, matching the registered fallback rule.
    by_title=defaultdict(list)
    for index,record in enumerate(records):
        title=record.get('normalized_title','')
        if title: by_title[title].append(index)
    for title,indices in by_title.items():
        components=defaultdict(list)
        for index in indices: components[find(index)].append(index)
        identified=[]; unidentified=[]
        for component in components.values():
            members=[records[i] for i in component]
            has_identifier=any(r.get('pmid') or r.get('doi') for r in members)
            (identified if has_identifier else unidentified).append(component)
        if len(identified)<=1:
            mergeable=identified+unidentified
            if mergeable:
                anchor=mergeable[0][0]
                for component in mergeable[1:]: union(anchor,component[0])
        else:
            # Multiple distinct PMID/DOI components make title-only matching
            # ambiguous. Retain them separately; merge only records with no IDs.
            if unidentified:
                anchor=unidentified[0][0]
                for component in unidentified[1:]: union(anchor,component[0])
    groups=defaultdict(list)
    for index,record in enumerate(records): groups[find(index)].append(record)
    out=[]
    priority={'PMID':0,'DOI':1,'TITLE':2,'RECORD':3}
    for members in groups.values():
        counts=defaultdict(int)
        for member in members:
            for key in dedup_keys(member):
                if key[0]!='TITLE': counts[key]+=1
        shared=[key for key,count in counts.items() if count>1]
        if not shared and len(members)>1 and len({m.get('pmid') or m.get('doi') for m in members})<=1:
            shared=[('TITLE',members[0].get('normalized_title',''))]
        rule=min(shared,key=lambda key:priority[key[0]])[0] if shared else dedup_keys(members[0])[0][0]
        out.append((rule,members))
    return out

def write_tsv(path: Path, fields, rows):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore'); w.writeheader(); w.writerows(rows)

def guard_review_outputs(out: Path):
    """Refuse regeneration once any reviewer or extraction data has been entered."""
    protected = [
        out/'screening/title_abstract_queue.tsv',
        out/'screening/fulltext_queue.tsv',
        out/'screening/fulltext_exclusions.tsv',
        out/'included_studies.tsv',
    ]
    stable_title_fields = {
        'screening_id','pmid','doi','year','title','abstract','query_id','duplicate_count','duplicate_sources'
    }
    for path in protected:
        if not path.is_file():
            continue
        try:
            with path.open(encoding='utf-8', newline='') as handle:
                reader = csv.DictReader(handle, delimiter='\t')
                if not reader.fieldnames:
                    raise SystemExit(f'Refusing rebuild: malformed review output {path}')
                stable = stable_title_fields if path.name == 'title_abstract_queue.tsv' else set()
                for row in reader:
                    if any((value or '').strip() for key, value in row.items() if key not in stable):
                        raise SystemExit(
                            'Refusing to overwrite reviewer/extraction data in '
                            f'{path.relative_to(out)}; archive it or implement an explicit merge first.'
                        )
        except (OSError, csv.Error) as exc:
            raise SystemExit(f'Refusing rebuild: cannot safely inspect {path}: {exc}') from exc

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--pubmed-dir',default='frailty_paper/review/pubmed'); ap.add_argument('--outdir',default='frailty_paper/review'); a=ap.parse_args()
    requested_base=Path(a.pubmed_dir)
    base=requested_base.resolve(); source_prefix=Path(requested_base.name)
    out=Path(a.outdir).resolve()
    guard_review_outputs(out)
    xmls=pubmed_xml_files(base)
    if not xmls: raise SystemExit(f'No PubMed XML files found under {base}')
    records=[]; file_rows=[]; parse_errors=[]
    for path in xmls:
        n=0
        try:
            for _, el in ET.iterparse(path,events=('end',)):
                if el.tag in {'PubmedArticle','PubmedBookArticle'}:
                    n+=1; records.append(parse_article(el,path,base,n,source_prefix)); el.clear()
        except Exception as e: parse_errors.append({'file':str(path),'error':str(e)})
        file_rows.append({'source_file':pubmed_source_path(path,base,source_prefix),'bytes':path.stat().st_size,'sha256':sha(path),'record_count':n})
    if parse_errors: raise SystemExit(f'XML parse failures: {parse_errors[:3]}')
    xml_record_count=len(records)
    manual_records_path=out/'manual_exports/manual_normalized_records.tsv'
    if manual_records_path.is_file():
        with manual_records_path.open(encoding='utf-8',newline='') as handle:
            reader=csv.DictReader(handle,delimiter='\t')
            required={'record_id','source_database','source_record_id','title','normalized_title','pmid','doi'}
            if not reader.fieldnames or not required.issubset(reader.fieldnames):
                raise SystemExit(f'Malformed normalized manual-export table: {manual_records_path}')
            for manual in reader:
                manual.setdefault('record_type','manual_database_record')
                manual.setdefault('query_id','manual_unknown')
                manual.setdefault('source_file','')
                manual.setdefault('abstract','')
                manual.setdefault('year','')
                manual.setdefault('journal','')
                manual.setdefault('authors','')
                manual.setdefault('import_status','READY' if manual.get('title') else 'MISSING_TITLE')
                records.append(manual)
    priority={'primary_sleep_frailty':0,'genetic_sleep_frailty':1,'frailty_neighborhood':2}
    dedup=[]; group_id=0
    duplicate_groups=group_duplicate_records(records)
    for rule, members in sorted(duplicate_groups,key=lambda group:min(dedup_keys(member) for member in group[1])):
        group_id+=1
        members.sort(key=lambda r:(priority.get(r['query_id'],9),r['source_file'],r['record_id']))
        retained=members[0]['record_id']
        for r in members:
            r['duplicate_group']=f'D{group_id:07d}'; r['retained_record_id']=retained; r['dedup_rule']=rule
        keep=dict(members[0]); keep['duplicate_count']=len(members); keep['duplicate_sources']=';'.join(sorted({m['source_file'] for m in members}))
        keep['duplicate_sources']=';'.join(sorted({f"{m.get('source_database','UNKNOWN')}:{m.get('source_record_id','')} ({m['source_file']})" for m in members}))
        dedup.append(keep)
    records.sort(key=lambda r:(priority.get(r['query_id'],9),r['source_file'],r['record_id']))
    all_fields=FIELDS
    write_tsv(out/'all_records.tsv',all_fields,records)
    dedup_fields=FIELDS+['duplicate_count','duplicate_sources']
    write_tsv(out/'deduplicated_records.tsv',dedup_fields,dedup)
    queue_fields=['screening_id','pmid','doi','year','title','abstract','query_id','duplicate_count','duplicate_sources','title_abstract_decision_r1','title_abstract_reason_r1','reviewer_1','title_abstract_decision_r2','title_abstract_reason_r2','reviewer_2','adjudication','adjudicator','notes']
    queue=[{'screening_id':r['record_id'],'pmid':r['pmid'],'doi':r['doi'],'year':r['year'],'title':r['title'],'abstract':r['abstract'],'query_id':r['query_id'],'duplicate_count':r['duplicate_count'],'duplicate_sources':r['duplicate_sources'],'title_abstract_decision_r1':'','title_abstract_reason_r1':'','reviewer_1':'','title_abstract_decision_r2':'','title_abstract_reason_r2':'','reviewer_2':'','adjudication':'','adjudicator':'','notes':''} for r in dedup]
    write_tsv(out/'screening/title_abstract_queue.tsv',queue_fields,queue)
    full_fields=['screening_id','pmid','doi','title','fulltext_location','fulltext_decision_r1','fulltext_reason_r1','reviewer_1','fulltext_decision_r2','fulltext_reason_r2','reviewer_2','adjudication','adjudicator','notes']
    write_tsv(out/'screening/fulltext_queue.tsv',full_fields,[])
    write_tsv(out/'screening/fulltext_exclusions.tsv',['screening_id','pmid','doi','title','primary_exclusion_reason','adjudication','reviewer','notes'],[])
    write_tsv(out/'included_studies.tsv',EXTRACTION,[])
    write_tsv(out/'pubmed_source_files.tsv',['source_file','bytes','sha256','record_count'],file_rows)
    query_counts=defaultdict(lambda:{'records':0,'unique_records':0})
    for r in records:
        if r.get('source_database')=='PubMed': query_counts[r['query_id']]['records']+=1
    for r in dedup:
        if r.get('source_database')=='PubMed': query_counts[r['query_id']]['unique_records']+=1
    count_rows=[{'query_id':q,'records_retrieved':v['records'],'unique_records_retained_by_source_priority':v['unique_records']} for q,v in sorted(query_counts.items())]
    write_tsv(out/'pubmed_acquisition_counts.tsv',['query_id','records_retrieved','unique_records_retained_by_source_priority'],count_rows)
    manual_count=sum(r.get('source_database')!='PubMed' for r in records)
    source_counts=defaultdict(lambda:{'records':0,'unique_records':0})
    for r in records: source_counts[r.get('source_database','UNKNOWN')]['records']+=1
    for r in dedup: source_counts[r.get('source_database','UNKNOWN')]['unique_records']+=1
    write_tsv(out/'source_record_counts.tsv',['database','records_imported_or_retrieved','unique_records_retained_by_priority'],[{'database':db,'records_imported_or_retrieved':v['records'],'unique_records_retained_by_priority':v['unique_records']} for db,v in sorted(source_counts.items())])
    prisma=[{'stage':'Records identified from PubMed searches','count':xml_record_count,'note':'Records in retained PubMed XML files'}, {'stage':'Records imported from other databases','count':manual_count,'note':'Untouched exports represented in the normalized manual-import table'}, {'stage':'Duplicate records removed across all sources','count':len(records)-len(dedup),'note':'PMID, then DOI, then normalized-title audit; every source record maps to a retained record'}, {'stage':'Records after deduplication','count':len(dedup),'note':'Unique records available for title/abstract screening'}, {'stage':'Records screened','count':0,'note':'Screening not started'}, {'stage':'Records excluded at title/abstract','count':0,'note':'Screening not started'}, {'stage':'Reports sought for retrieval','count':0,'note':'Screening not started'}, {'stage':'Reports not retrieved','count':0,'note':'Screening not started'}, {'stage':'Reports assessed for eligibility','count':0,'note':'Screening not started'}, {'stage':'Reports excluded after full text','count':0,'note':'Screening not started'}, {'stage':'Studies included','count':0,'note':'No studies included yet'}]
    write_tsv(out/'prisma_flow.tsv',['stage','count','note'],prisma)
    output_names=['all_records.tsv','deduplicated_records.tsv','screening/title_abstract_queue.tsv','screening/fulltext_queue.tsv','screening/fulltext_exclusions.tsv','included_studies.tsv','pubmed_source_files.tsv','pubmed_acquisition_counts.tsv','source_record_counts.tsv','prisma_flow.tsv']
    output_rows=[]
    for name in output_names:
        p=out/name
        output_rows.append({'file':str(p.relative_to(out)),'bytes':p.stat().st_size,'sha256':sha(p)})
    write_tsv(out/'record_build_outputs.tsv',['file','bytes','sha256'],output_rows)
    manifest={'built_at_utc':datetime.now(timezone.utc).isoformat(),'command':'python frailty_paper/scripts/08_build_review_records.py','script_sha256':sha(Path(__file__)),'xml_files':len(xmls),'xml_records':xml_record_count,'manual_import_records':manual_count,'total_source_records':len(records),'record_type_counts':{'journal_article':sum(r['record_type']=='journal_article' for r in records),'book_article':sum(r['record_type']=='book_article' for r in records),'manual_database_record':sum(r['record_type']=='manual_database_record' for r in records)},'unique_records':len(dedup),'duplicate_occurrences':len(records)-len(dedup),'dedup_rule':'Exact PMID; normalized DOI where titles agree or a unique title candidate exists; normalized title only without conflicting IDs; otherwise unique source record','retained_source_priority':['primary_sleep_frailty','genetic_sleep_frailty','frailty_neighborhood','manual_database_sources'],'source_manifest':'pubmed_source_files.tsv','output_manifest':'record_build_outputs.tsv','output_manifest_sha256':sha(out/'record_build_outputs.tsv')}
    (out/'record_build_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(manifest,indent=2))
if __name__=='__main__': main()
