"""Pinned author training identities, sequences and folds; omit activity labels."""
import urllib.request,ssl,hashlib,gzip,csv,json,concurrent.futures
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'outputs/sleep_frailty_ai_v9_validation/next_source_screen/candidate_pretrained_readers'
COMMIT='a20dbc7d0b08daae3af1aebcfb8a71d9f198ab5d'


def main():
    def fetch(cell):
        path=OUT/f'{cell}_original_pretraining_identity_and_folds.tsv.gz';record=OUT/f'{cell}_original_pretraining_manifest_receipt.json'
        if record.exists():
            r=json.loads(record.read_text());assert hashlib.sha256(path.read_bytes()).hexdigest()==r['manifest_sha256'];return r
        url=f'https://raw.githubusercontent.com/autosome-ru/human_legnet/{COMMIT}/datasets/original/{cell}.tsv';raw_hash=hashlib.sha256();nbytes=0;nrows=0
        if cell=='WTC11':handle=Path('/Volumes/Extreme SSD/codex-sleep-frailty-work/sleep_frailty_ml_v8/sources/legnet_datasets_original_WTC11.tsv').open('rb')
        else:handle=urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'SleepFrailtyResearch/1.0'}),context=ssl.create_default_context(cafile='/etc/ssl/cert.pem'),timeout=45)
        with handle,gzip.open(path,'wt') as output:
            header_line=handle.readline();raw_hash.update(header_line);nbytes+=len(header_line);header=header_line.decode().rstrip('\r\n').split('\t');assert set(header)=={'seq_id','seq','mean_value','fold','rev'}
            keep=['seq_id','seq','fold','rev'];columns=[header.index(k) for k in keep];writer=csv.writer(output,delimiter='\t');writer.writerow(keep)
            for line in handle:
                nbytes+=len(line);assert nbytes<=200_000_000;raw_hash.update(line);fields=line.decode().rstrip('\r\n').split('\t');assert len(fields)==len(header)
                row=[fields[i] for i in columns];assert len(row[1])==230 and not set(row[1])-set('ACGT') and int(row[2]) in range(1,11) and row[3] in ['0','1'];writer.writerow(row);nrows+=1
        r={'cell':cell,'source_url':url,'repository_commit':COMMIT,'original_input_bytes':nbytes,'original_input_sha256':raw_hash.hexdigest(),'rows':nrows,'manifest_path':str(path),'manifest_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
           'columns':['seq_id','seq','fold','rev'],'activity_labels_omitted_from_saved_manifest':True,'new_experimental_effects_or_model_performance_used':False}
        record.write_text(json.dumps(r,indent=2)+'\n');print('AUTHOR PRETRAINING FOLDS VERIFIED',cell,nrows,flush=True);return r
    with concurrent.futures.ThreadPoolExecutor(2) as pool:records=list(pool.map(fetch,['HepG2','K562','WTC11']))
    (OUT/'all_pretraining_fold_manifest_receipt.json').write_text(json.dumps({'records':records,'purpose':'Check actual reference-activity pretraining exposure; no new-study outcomes or activity-based row selection.'},indent=2)+'\n')


if __name__=='__main__':main()
