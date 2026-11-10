"""Original released count tables for auditing already revealed training targets."""
import zipfile
from channel_common import *
from fetch_cell_matched_pretrained_models import RemoteArchive

def main():
    assert_frozen_previous();folder=DEST/'source_candidates';folder.mkdir(exist_ok=True);url='https://zenodo.org/api/records/15297965/files/data_MPRA.zip/content';remote=RemoteArchive(url,929286104);wanted=['OL42_20221220_counts.out.gz','OL42_20221220_normalized_counts.out.gz'];records=[]
    with zipfile.ZipFile(remote) as archive:
        for name in wanted:
            member='data/MPRA/counts/'+name;info=archive.getinfo(member);assert info.file_size<45_000_000;path=folder/name
            if not path.exists():
                partial=path.with_suffix('.partial');n=0
                try:
                    with archive.open(info) as source,partial.open('wb') as handle:
                        while block:=source.read(1_000_000):handle.write(block);n+=len(block)
                    assert n==info.file_size;partial.replace(path)
                except BaseException:
                    if partial.exists():partial.unlink()
                    raise
            # Separate gzip CRC verifies the released inner compressed table.
            with gzip.open(path,'rb') as handle:
                while handle.read(1_000_000):pass
            records.append({'member':member,'path':str(path),'bytes':path.stat().st_size,'sha256':sha(path),'zip_and_gzip_CRC_verified':True});print('ORIGINAL TRAINING COUNTS RETRIEVED',name,flush=True)
            with gzip.open(path,'rt') as handle:print('COUNT TABLE HEADER',repr(handle.readline()[:3000]),repr(handle.readline()[:900]),flush=True)
    write_json(folder/'Siraj_original_count_retrieval_receipt.json',{'source':url,'files':records,'range_bytes_read':remote.total,'role':'Audit already revealed source identities, A/B mappings and replicate noise. Not a new independent model validation or an inferred original GLM sign convention.'});assert_frozen_previous()

if __name__=='__main__':main()
