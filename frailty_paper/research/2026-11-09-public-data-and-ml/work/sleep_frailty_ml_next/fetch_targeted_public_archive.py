"""Bounded HTTP-range retrieval of small, named members of a public ZIP64 archive."""
import io,json,hashlib,urllib.request,ssl,zipfile
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'outputs/sleep_frailty_ai_v9_validation/next_source_screen'


class BoundedRemoteZip(io.RawIOBase):
    def __init__(self,url,size):self.url=url;self.size=size;self.position=0;self.total=0;self.ctx=ssl.create_default_context(cafile='/etc/ssl/cert.pem')
    def seekable(self):return True
    def readable(self):return True
    def tell(self):return self.position
    def seek(self,offset,whence=0):
        self.position=offset if whence==0 else self.position+offset if whence==1 else self.size+offset
        assert 0<=self.position<=self.size;return self.position
    def read(self,size=-1):
        if size<0:size=self.size-self.position
        size=min(size,self.size-self.position)
        if size==0:return b''
        assert size<=2_000_000 and self.total+size<=10_000_000,'Metadata/member transfer bound exceeded'
        request=urllib.request.Request(self.url,headers={'Range':f'bytes={self.position}-{self.position+size-1}','User-Agent':'SleepFrailtyResearch/1.0'})
        with urllib.request.urlopen(request,context=self.ctx,timeout=45) as response:
            assert response.status==206 and response.headers.get('Content-Range','').startswith(f'bytes {self.position}-'),'Server did not honor bounded range request'
            data=response.read(size+1)
        assert len(data)==size;self.position+=size;self.total+=size;return data


def main():
    metadata=json.loads((OUT/'Tardaguila2026_Zenodo_metadata.json').read_text());archive=metadata['files'][0];remote=BoundedRemoteZip(archive['links']['self'],archive['size']);base='Zenodo_V2F_paper/'
    wanted={
        'README.md':'Tardaguila2026_README.txt',
        '04_data_analysis/MPRA/Analysis/README.md':'Tardaguila2026_MPRA_analysis_README.txt',
        '04_data_analysis/MPRA/Alignment_and_counts/README.md':'Tardaguila2026_MPRA_alignment_README.txt',
        '04_data_analysis/MPRA/Alignment_and_counts/Dependencies/Library_TRIMMED_15bp_Carried_Variants.fasta':'Tardaguila2026_actual_library.fasta',
        '04_data_analysis/MPRA/Analysis/Dependencies/MPRA_sequence_information.tsv':'Tardaguila2026_sequence_information.tsv',
        '01_Tables/Table_S3/Table_S3.tsv':'Tardaguila2026_TableS3.tsv',
        '01_Tables/Table_S4/Table_S4.tsv':'Tardaguila2026_TableS4.tsv',
        '01_Tables/Table_S7/Table_S7_Oligos_and_sgRNAs_screenings_v2.xlsx':'Tardaguila2026_TableS7.xlsx'}
    receipt=[]
    with zipfile.ZipFile(remote) as z:
        for member,name in wanted.items():
            info=z.getinfo(base+member);assert info.file_size<=2_000_000 and info.compress_size<=1_000_000
            content=z.read(info);assert len(content)==info.file_size # zipfile verifies the stored CRC
            path=OUT/name
            if path.exists():assert path.read_bytes()==content
            else:path.write_bytes(content)
            receipt.append({'member':info.filename,'path':str(path),'bytes':len(content),'sha256':hashlib.sha256(content).hexdigest(),'archive_CRC_verified':True})
    (OUT/'targeted_archive_retrieval_receipt.json').write_text(json.dumps({'source_doi':'10.5281/zenodo.17860579','article_doi':'10.1016/j.celrep.2026.116937','archive_bytes':archive['size'],'actual_network_range_bytes':remote.total,'files':receipt,'new_model_training_or_performance_evaluation':False},indent=2)+'\n')
    print('TARGETED FREE SOURCE FILES VERIFIED',len(receipt),'NETWORK BYTES',remote.total,flush=True)


if __name__=='__main__':main()
