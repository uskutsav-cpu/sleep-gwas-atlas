"""Small selected author checkpoints, retrieved by bounded ranges, not deserialized."""
import io,json,hashlib,urllib.request,ssl,zipfile
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'outputs/sleep_frailty_ai_v9_validation/next_source_screen/candidate_pretrained_readers'


class RemoteArchive(io.RawIOBase):
    def __init__(self,url,size):self.url=url;self.size=size;self.position=0;self.total=0;self.ctx=ssl.create_default_context(cafile='/etc/ssl/cert.pem')
    def readable(self):return True
    def seekable(self):return True
    def tell(self):return self.position
    def seek(self,offset,whence=0):
        self.position=offset if whence==0 else self.position+offset if whence==1 else self.size+offset
        assert 0<=self.position<=self.size;return self.position
    def read(self,size=-1):
        if size<0:size=self.size-self.position
        size=min(size,self.size-self.position)
        assert size<=20_000_000 and self.total+size<=60_000_000
        chunks=[]
        while size:
            n=min(size,2_000_000);request=urllib.request.Request(self.url,headers={'Range':f'bytes={self.position}-{self.position+n-1}','User-Agent':'SleepFrailtyResearch/1.0'})
            with urllib.request.urlopen(request,context=self.ctx,timeout=45) as response:
                assert response.status==206 and response.headers.get('Content-Range','').startswith(f'bytes {self.position}-');chunk=response.read(n+1)
            assert len(chunk)==n;chunks.append(chunk);self.position+=n;self.total+=n;size-=n
        return b''.join(chunks)


def main():
    metadata=json.loads((OUT/'Zenodo_10558183_metadata.json').read_text());f=next(r for r in metadata['files'] if r['key']=='final_dump.zip');remote=RemoteArchive(f['links']['self'],f['size']);receipt=[]
    with zipfile.ZipFile(remote) as z:
        for cell in ['HepG2','K562','WTC11']:
            base=f'final_dump/models/{cell}/md_shift_reverse_noavg_noch/'
            for member,target in [('config.json',cell+'_author_config.json'),('best_model_test1_val2.ckpt',cell+'_author_test1_val2.ckpt')]:
                info=z.getinfo(base+member);assert info.file_size<=20_000_000;content=z.read(info);path=OUT/target
                if path.exists():assert path.read_bytes()==content
                else:path.write_bytes(content)
                receipt.append({'source_member':info.filename,'path':str(path),'bytes':len(content),'sha256':hashlib.sha256(content).hexdigest(),'archive_CRC_verified':True})
            print('AUTHOR CHECKPOINT RETRIEVED',cell,flush=True)
    (OUT/'cell_matched_pretrained_retrieval_receipt.json').write_text(json.dumps({'source_doi':'10.5281/zenodo.10558183','article_doi':'10.1038/s41586-024-08430-9','fixed_selection':'test1/validation2 for all three cells, chosen by matching filename without inspecting any new outcome or model performance. This is source access and architecture verification, not the paper-recommended full 90-model ensemble.',
        'range_bytes_transferred':remote.total,'files':receipt,'checkpoints_deserialized':False,'new_benchmark_performance_calculated':False},indent=2)+'\n')
    print('SELECTED AUTHOR MODELS RETRIEVED; TENSOR-SAFETY AND PRETRAINING-OVERLAP AUDITS STILL REQUIRED',flush=True)


if __name__=='__main__':main()
