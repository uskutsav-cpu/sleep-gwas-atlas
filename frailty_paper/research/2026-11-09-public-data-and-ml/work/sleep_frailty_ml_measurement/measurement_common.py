"""Measurement-aware development, preserving all prior failed experiments."""
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'work/sleep_frailty_ml_channels'))
from channel_common import ROOT,OLD,assert_frozen_previous,sha,write_json,write_csv,csv_rows,rc,BASE,V6,V7,V8,DATA,np,json,collections,csv,gzip,hashlib,os
CHANNEL=ROOT/'outputs/sleep_frailty_ml_channels'
SOURCE=CHANNEL/'source_candidates'
DEST=ROOT/'outputs/sleep_frailty_ml_measurement'
DEST.mkdir(exist_ok=True)

def assert_prior_versions():
    assert_frozen_previous()
    for name,h in json.loads((CHANNEL/'final_progress_verification.json').read_text())['files'].items():assert sha(CHANNEL/name)==h,name
