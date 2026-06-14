"""Spatially correlated fibrosis as a low-coupling conductivity field.

White noise smoothed by a Gaussian kernel (correlation length corr_len, in cells)
is thresholded at the quantile that yields the requested area fraction. Smoothing
gives realistic connected patches instead of isolated squares; thresholding pins
the density exactly. Returns a coupling-scale field (d_scar inside, d_healthy out).
"""
from __future__ import annotations

import numpy as np
from scipy.ndimage import gaussian_filter


def correlated_fibrosis(shape, density=0.0, corr_len=4.0, d_healthy=1.0,
                        d_scar=0.05, seed=0):
    field = np.full(shape, d_healthy, dtype=float)
    if density <= 0:
        return field
    rng = np.random.default_rng(seed)
    noise = rng.standard_normal(shape)
    smooth = gaussian_filter(noise, sigma=corr_len, mode="reflect")
    thresh = np.quantile(smooth, density)
    field[smooth <= thresh] = d_scar
    return field
