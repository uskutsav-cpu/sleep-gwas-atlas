from pathlib import Path
import json,csv,hashlib,numpy as np
ROOT=Path(__file__).resolve().parents[2];WORK=Path(__file__).resolve().parent;DATA=Path('/Volumes/Extreme SSD/codex-sleep-frailty-work/sleep_frailty_ml_v4');OUT=DATA/'outputs';V3=ROOT/'outputs/sleep_frailty_ai_v3';OLD=ROOT/'work/sleep_frailty_ml';INDEX=DATA.parent/'stage4_atlas_cache/avi_snvs.tsv.gz.tbi'
AVI_URL='https://grr.iossifovlab.com/hg38/scores/AVI/alphagenome_variant_impact_score_snvs.tsv.gz'
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()
def write_json(path,value):
    p=Path(path);temp=p.with_suffix(p.suffix+'.tmp');temp.write_text(json.dumps(value,indent=2)+'\n');temp.replace(p)
