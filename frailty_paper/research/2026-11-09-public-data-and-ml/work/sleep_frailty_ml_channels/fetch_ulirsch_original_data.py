"""Bounded retrieval from the original study's official GEO FTP data deposit."""
import urllib.request,time
from channel_common import *

def main():
    dest=DEST/'source_candidates';dest.mkdir(exist_ok=True);base='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE87nnn/GSE87711/';files=[('suppl/','GSE87711_RBC_MPRA_minP_raw.txt.gz'),('suppl/','GSE87711_design.txt.gz'),('suppl/','GSE87711_design_barcodes.txt.gz'),('soft/','GSE87711_family.soft.gz')];records=[]
    for folder,name in files:
        url=base+folder+name;path=dest/name
        if not path.exists():
            with urllib.request.urlopen(url,timeout=45) as response:data=response.read(8_000_001)
            assert len(data)<=8_000_000
            # Full gzip decompression verifies CRC; data are not executed.
            with gzip.GzipFile(fileobj=__import__('io').BytesIO(data)) as handle:
                while handle.read(1_000_000):pass
            path.write_bytes(data)
        records.append({'url':url,'path':str(path),'bytes':path.stat().st_size,'sha256':sha(path),'gzip_CRC_verified':True});print('ORIGINAL ULIRSCH FILE RETRIEVED',name,path.stat().st_size,flush=True)
    write_json(dest/'Ulirsch_original_retrieval_receipt.json',{'source':'Ulirsch et al. 2016, doi:10.1016/j.cell.2016.04.048; original GEO GSE87711','files':records,
              'sequence_or_count_files_not_model_generated':True,'full_raw_FASTQ_or_large_archive_downloaded':False,'experimental_effects_fitted':False,'model_predictions_or_performance_calculated':False,'candidate_only_not_yet_eligible_validation':True})
    for folder,name in files[:3]:
        with gzip.open(dest/name,'rt') as handle:print('SOURCE HEADER',name,repr(handle.readline()[:1800]),repr(handle.readline()[:800]),flush=True)

if __name__=='__main__':main()
