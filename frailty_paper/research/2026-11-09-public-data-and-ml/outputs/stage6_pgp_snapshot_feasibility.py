#!/usr/bin/env python3
"""Aggregate feasibility only; do not emit participant identifiers or records."""
from pathlib import Path
import sqlite3,gzip,hashlib,json,requests
W=Path('work');O=Path('outputs');p=W/'stage6_pgp_untap_snapshot.sqlite3.gz';raw=gzip.open(p,'rb').read();c=sqlite3.connect(':memory:');c.deserialize(raw)
tables={}
for t in ['demographics','survey','conditions','test_results','uploaded_data']:
 tables[t]={'rows':c.execute('select count(*) from '+t).fetchone()[0],'distinct_people':c.execute('select count(distinct human_id) from '+t).fetchone()[0]}
q='select count(distinct a.human_id) from conditions a join uploaded_data b on a.human_id=b.human_id where lower(a.name) like ? and b.data_type=? and b.download_url is not null and length(trim(b.download_url))>0';typed={}
for typ in ['23andMe','Complete Genomics','Illumina']:
 typed[typ]={'listed_url_people':c.execute('select count(distinct human_id) from uploaded_data where data_type=? and download_url is not null and length(trim(download_url))>0',(typ,)).fetchone()[0],'insomnia_label_people':c.execute(q,('%insom%',typ)).fetchone()[0]}
union=c.execute("select count(distinct a.human_id) from conditions a join uploaded_data b on a.human_id=b.human_id join demographics d on a.human_id=d.human_id where lower(a.name) like '%insom%' and b.data_type in ('23andMe','Complete Genomics','Illumina') and b.download_url is not null and length(trim(b.download_url))>0").fetchone()[0]
url='https://download.ac2it.arvadosapi.com/c=3ccedff14c04af4d96e41e96648c55e2+3173/hu-pgp.sqlite3.gz';res={'url':url,'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'publisher_checksum_verified':False,'uncompressed_bytes':len(raw),'snapshot_date':'not established','tables':tables,'types_not_mutually_exclusive':typed,'insomnia_labeled_with_demographics_and_any_listed_download':union,'caution':'Condition name substring is not a validated diagnosis, missing condition entries are not negative controls, listed URL is not verified usable genotype data, participant overlap and availability of score alleles and frailty outcomes remain unknown. Historical public snapshot may not reflect current participant sharing choices.'};(O/'stage6_pgp_snapshot_feasibility.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res,indent=2))
