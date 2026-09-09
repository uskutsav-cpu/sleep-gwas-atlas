from pathlib import Path
import pytest

@pytest.mark.parametrize('name',['common.R','placo.R','susie.R','coloc.R','lava.R','genomicsem.R','mr.R','check_native.R','native_environment.R','native_cli.py'])
def test_bundled_native_adapter_matches_source(name):
    root=Path(__file__).resolve().parents[1]
    assert (root/'scripts'/name).read_bytes()==(root/'brain6/resources'/name).read_bytes()
