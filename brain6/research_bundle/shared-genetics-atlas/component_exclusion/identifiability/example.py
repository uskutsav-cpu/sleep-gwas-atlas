#!/usr/bin/env python3
"""Two hypothetical correlation matrices consistent with one composite result."""
import math
observed=0.395
for sleep_exhaustion in (0.2,0.6):
    # P and E have equal genetic variance, are standardized and uncorrelated.
    sleep_physical=observed*math.sqrt(2)-sleep_exhaustion
    determinant=1-sleep_physical**2-sleep_exhaustion**2
    reconstructed=(sleep_physical+sleep_exhaustion)/math.sqrt(2)
    assert determinant>0 and abs(reconstructed-observed)<1e-12
    print(f'r(S,E)={sleep_exhaustion:.4f}, r(S,P)={sleep_physical:.4f}, r(S,F)={reconstructed:.4f}, determinant={determinant:.4f}')
