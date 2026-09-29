import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "brain6/scripts"))

from figure_style import best_contrast_text
from build_locked_global import derive_sleep_traits
import pandas as pd
import pytest


def test_annotation_text_uses_accessible_contrast_on_dark_and_light_cells():
    assert best_contrast_text("#000000") == "white"
    assert best_contrast_text("#ffffff") == "black"


def test_magma_heatmap_extremes_receive_readable_text():
    from matplotlib import colormaps

    magma = colormaps["magma"]
    assert best_contrast_text(magma(0.0)) == "white"
    assert best_contrast_text(magma(1.0)) == "black"


def test_sleep_traits_are_derived_from_panel_without_a_fixed_count():
    panel = pd.DataFrame({
        "domain": ["sleep", "other", "sleep", "sleep"],
        "trait_id": ["sleep_a", "disorder", "sleep_b", "sleep_c"],
    })
    assert derive_sleep_traits(panel) == ["sleep_a", "sleep_b", "sleep_c"]


def test_sleep_panel_rejects_duplicate_or_missing_trait_ids():
    duplicate = pd.DataFrame({"domain": ["sleep", "sleep"], "trait_id": ["insomnia", "insomnia"]})
    with pytest.raises(ValueError, match="unique"):
        derive_sleep_traits(duplicate)
    empty = pd.DataFrame({"domain": ["other"], "trait_id": ["mdd"]})
    with pytest.raises(ValueError, match="unique, nonempty"):
        derive_sleep_traits(empty)
