"""Fetch small public author source tables; verify exact Git blob bytes."""
import urllib.request,ssl,certifi,time
from representation_common import *


def blob_sha(data):return hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()


def main():
    assert_versions_all();folder=DEST/'neural_source_trace';m=json.loads((folder/'source_metadata.json').read_text());rows=[]
    for item in m['files']:
        path=folder/Path(item['path']).name
        if path.exists():
            raw=path.read_bytes()
            if blob_sha(raw)!=item['sha'] and blob_sha(raw[:-1])==item['sha']:
                raw=raw[:-1];path.write_bytes(raw)
        else:
            url=f'https://raw.githubusercontent.com/{m["repository"]}/{m["commit"]}/{item["path"]}'
            with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'SleepFrailtyResearch/1.0'}),context=ssl.create_default_context(cafile=certifi.where()),timeout=60) as response:
                raw=response.read()
            assert len(raw)==item['size'] and blob_sha(raw)==item['sha'];path.write_bytes(raw)
        assert len(raw)==item['size'] and blob_sha(raw)==item['sha']
        if path.suffix=='.gz':
            with gzip.open(path,'rb') as handle:
                while handle.read(1024*1024):pass
        rows.append({'path':str(path),'bytes':len(raw),'git_blob_sha1':blob_sha(raw),'sha256':sha(path)})
        print('ORIGINAL NEURAL TABLE VERIFIED',path.name,len(raw),flush=True)
    write_json(folder/'retrieval_receipt.json',{'passed':True,'repository':m['repository'],'commit':m['commit'],'files':rows,
       'original_counts_or_GLM_refitted':False,'new_outside_test_read':False})
    assert_versions_all()


if __name__=='__main__':main()
