"""Regression checks against specific scientific overclaims."""
import importlib.util
from pathlib import Path
import pytest
p=Path(__file__).resolve().parents[1]/'scripts/validate_manuscript.py'
import sys
sys.path.insert(0,str(p.parent))
spec=importlib.util.spec_from_file_location('claimcheck',p);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
TEXT=(p.parent.parent/'MANUSCRIPT/brain6_manuscript.tex').read_text()
def test_current_required_limits():assert m.validate_text(TEXT)==[]
@pytest.mark.parametrize('claim',['we identified a novel shared causal variant','we established independent two-trait replication','canonical LAVA passed QC','all LD matrices passed','a causal effector gene was confirmed'])
def test_rejects_headline_overclaim(claim):assert any('UNSUPPORTED_POSITIVE_CLAIM' in x for x in m.validate_text(TEXT+'\n'+claim))
def test_failed_qc_cannot_disappear():assert any('MISSING_REQUIRED' in x for x in m.validate_text(TEXT.replace('FAILED\\_QC\\_NOT\\_PROMOTED','PROMOTED')))
