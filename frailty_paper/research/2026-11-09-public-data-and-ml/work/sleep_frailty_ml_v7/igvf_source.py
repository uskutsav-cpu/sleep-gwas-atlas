"""Public released IGVF files through the documented API and publisher checksums."""
from common import *
def metadata(accession):
    p=SOURCES/(accession+'_metadata.json')
    if not p.exists():fetch('https://api.data.igvf.org/'+accession+'/?format=json',p,8_000_000)
    d=json.loads(p.read_text());assert d['status']=='released',accession
    return d
def download(accession,limit=60_000_000):
    requested=accession;chain=[]
    while True:
        p=SOURCES/(accession+'_metadata.json')
        if not p.exists():fetch('https://api.data.igvf.org/'+accession+'/?format=json',p,8_000_000)
        d=json.loads(p.read_text());chain.append({'accession':accession,'status':d['status'],'metadata_sha256':sha(p)})
        if d['status']=='released':break
        targets=d.get('superseded_by',[])
        assert d['status'] in ['revoked','archived','superseded'] and len(targets)==1,'No unique released replacement for '+accession
        accession=targets[0].strip('/').split('/')[-1];assert accession not in [r['accession'] for r in chain]
    assert not d.get('controlled_access',False)
    p=SOURCES/(accession+'.gz')
    uri=d.get('s3_uri','')
    # Newer released files use the API's public download view. Do not guess
    # private-bucket URLs, require an account, or bypass access controls.
    if uri.startswith('s3://igvf-public/'):url='https://igvf-public.s3.amazonaws.com/'+uri.split('igvf-public/')[1]
    else:url='https://api.data.igvf.org'+d['href']
    if not p.exists():fetch(url,p,limit,b'\x1f\x8b')
    h=hashlib.md5()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    assert h.hexdigest()==d['md5sum'],accession
    receipt={'requested_accession':requested,'accession':accession,'replacement_chain':chain,'url':url,'path':str(p),'metadata_sha256':sha(SOURCES/(accession+'_metadata.json')),'bytes':p.stat().st_size,'publisher_MD5_matches':True,'sha256':sha(p),'status':d['status'],'content_type':d.get('content_type'),'submitted_file_name':d.get('submitted_file_name')}
    write_json(OUT/(accession+'_source_receipt.json'),receipt)
    if requested!=accession:write_json(OUT/(requested+'_replacement_receipt.json'),receipt)
    return p
if __name__=='__main__':
    import sys
    for acc in sys.argv[1:]:print(download(acc),flush=True)
