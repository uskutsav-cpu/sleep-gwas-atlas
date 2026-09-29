"""Small styling helpers shared by Brain6 publication figures and tests."""
from __future__ import annotations

import numpy as np
import matplotlib


def best_contrast_text(background) -> str:
    """Choose black or white text by WCAG contrast against an RGB/RGBA color."""
    srgb = np.asarray(matplotlib.colors.to_rgb(background), dtype=float)
    linear = np.where(srgb <= 0.04045, srgb / 12.92, ((srgb + 0.055) / 1.055) ** 2.4)
    luminance = float(np.dot(linear, [0.2126, 0.7152, 0.0722]))
    white_contrast = 1.05 / (luminance + 0.05)
    black_contrast = (luminance + 0.05) / 0.05
    return "white" if white_contrast >= black_contrast else "black"
