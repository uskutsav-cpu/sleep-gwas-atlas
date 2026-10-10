#!/usr/bin/env python3
"""Bounded HEAD-only original-release metadata; never read response bodies."""
import csv
import datetime
import hashlib
import json
from pathlib import Path
import ssl
import urllib.error
import urllib.parse
import urllib.request

P=Path(__file__).resolve().parents[1]
ROOT=P.parent


class HeadOnlyRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self,request,fp,code,message,headers,newurl):
        if request.get_method()!='HEAD' or urllib.parse.urlparse(newurl).scheme!='https':
            raise RuntimeError('NON_HEAD_OR_NON_HTTPS_REDIRECT_REFUSED')
        return urllib.request.Request(newurl,headers={'User-Agent':'sleep-validation-provenance-design/1.0'},method='HEAD')


def main():
    queue=ROOT/'discovery_extension/results/replication/replication_source_queue.tsv'
    with queue.open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
    unique={}
    for row in rows:
        if row['source_curation_status']=='COMPLETE_BEFORE_RESULTS':unique.setdefault(row['replication_source_id'],row)
    if len(unique)!=13:raise RuntimeError('ORIGINAL13_CARDINALITY_DIFFERS')
    import certifi
    context=ssl.create_default_context(cafile=certifi.where())
    opener=urllib.request.build_opener(HeadOnlyRedirect(),urllib.request.HTTPSHandler(context=context))
    receipts=[]
    for index,(source,row) in enumerate(unique.items(),1):
        r=dict(index=index,source_id=source,original_url=row['replication_source_url'],method='HEAD',timeout_seconds=10,
               recorded_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),GWAS_body_read_or_downloaded=False)
        try:
            request=urllib.request.Request(row['replication_source_url'],headers={'User-Agent':'sleep-validation-provenance-design/1.0'},method='HEAD')
            with opener.open(request,timeout=10) as response:
                r.update(status=response.status,final_url=response.geturl(),headers={k.lower():v for k,v in response.headers.items()})
                if len(json.dumps(r['headers']))>256*(1<<10):raise RuntimeError('METADATA_HEADERS_EXCEED256KiB')
            h=r['headers'];checks=dict(length=h.get('content-length')==row['replication_content_length_bytes'],etag=h.get('etag','').strip('"')==row['replication_etag'])
            if row['replication_storage_mode']=='VERSIONED_REMOTE_STREAMING':checks['generation']=h.get('x-goog-generation')==row['replication_source_generation']
            r['identity_checks']=checks;r['verdict']='HEAD_IDENTITY_METADATA_MATCH_ONLY' if all(checks.values()) else 'HEAD_METADATA_IDENTITY_DIFFERS_REQUIRES_REVIEW'
        except Exception as error:
            r.update(verdict='HEAD_METADATA_UNAVAILABLE_NO_BODY_FALLBACK',error=type(error).__name__+': '+str(error))
        receipts.append(r);print(json.dumps(dict(index=index,source_id=source,verdict=r['verdict'])),flush=True)
    out=Path(__file__).with_suffix('.json')
    record=dict(scope='Original13 release HEAD metadata only; no raw/processed GWAS read, decompression, acquisition, harmonization, munging, merge or fit.',queue_sha256=hashlib.sha256(queue.read_bytes()).hexdigest(),collector_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),receipts=receipts)
    with out.open('x') as f:json.dump(record,f,indent=2);f.write('\n')
    print(json.dumps(dict(receipt=str(out),sha256=hashlib.sha256(out.read_bytes()).hexdigest(),sources=len(receipts),GWAS_body_bytes_read=0)),flush=True)


if __name__=='__main__':main()
