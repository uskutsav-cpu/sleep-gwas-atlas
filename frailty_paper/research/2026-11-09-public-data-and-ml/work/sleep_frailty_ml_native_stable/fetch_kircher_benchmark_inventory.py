"""Original-author OSF identities for four published saturation benchmarks."""
import csv,hashlib,io,json,ssl,urllib.request
from pathlib import Path
import certifi

ROOT=Path(__file__).resolve().parents[2]
DEST=ROOT/'outputs/sleep_frailty_ml_native_stable/benchmark_exposure/kircher'
LISTING='https://api.osf.io/v2/nodes/75b2m/files/osfstorage/5bfc218d576a4f0017b38543/?page[size]=100'


def main():
    DEST.mkdir(exist_ok=True);ctx=ssl.create_default_context(cafile=certifi.where())
    with urllib.request.urlopen(LISTING,context=ctx,timeout=30) as response:raw=response.read()
    listing=json.loads(raw);assert listing['links']['next'] is None
    selected=[x for x in listing['data'] if any(x['attributes']['name'].startswith('GRCh38_'+gene) for gene in ['F9.','LDLR','PKLR','SORT1'])]
    assert len(selected)>=4
    (DEST/'OSF_GRCh38_original_file_listing.json').write_bytes(raw);records=[]
    for item in selected:
        a=item['attributes'];assert a['kind']=='file' and a['size']<250000
        url=item['links']['download']
        with urllib.request.urlopen(url,context=ctx,timeout=30) as response:content=response.read()
        h=hashlib.sha256(content).hexdigest();assert len(content)==a['size'] and h==a['extra']['hashes']['sha256']
        path=DEST/a['name'];path.write_bytes(content)
        lines=content.decode().splitlines();print(a['name'],len(content),'header',lines[0],flush=True)
        records.append({'name':a['name'],'url':url,'OSF_id':item['id'],'OSF_version':a['current_version'],
             'sha256':h,'bytes':len(content),'date_modified':a['date_modified'],'header':lines[0],'first_identity_row':lines[1].split('\t')[:5]})
    (DEST/'original_saturation_inventory_retrieval.json').write_text(json.dumps({'passed':True,'source_paper':'https://pmc.ncbi.nlm.nih.gov/articles/PMC6687891/',
       'reader_benchmark_paper':'https://pmc.ncbi.nlm.nih.gov/articles/PMC11903340/','OSF_listing_url':LISTING,
       'OSF_listing_sha256':hashlib.sha256(raw).hexdigest(),'files':records,'reserved_new_study_outcomes_read':False,
       'raw_file_coefficients_not_analyzed':True},indent=2)+'\n')


if __name__=='__main__':main()
