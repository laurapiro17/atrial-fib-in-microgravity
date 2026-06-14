import numpy as np
from afib_microgravity.metrics import is_sustained_af, ps_density, bootstrap_ci

def test_sustained_af_true_when_ps_persists():
    series = [5, 6, 4, 7, 5, 6]
    assert is_sustained_af(series, window_frac=0.5, threshold=2) is True

def test_sustained_af_false_when_quiescent():
    series = [3, 2, 1, 0, 0, 0]
    assert is_sustained_af(series, window_frac=0.5, threshold=2) is False

def test_ps_density_normalises_by_area():
    assert ps_density(8, area_cells=40000) == 8 / 40000 * 1e4

def test_bootstrap_ci_brackets_mean():
    vals = [1.0, 2.0, 3.0, 4.0, 5.0]
    mean, lo, hi = bootstrap_ci(vals, n_boot=500, seed=0)
    assert abs(mean - 3.0) < 1e-9
    assert lo <= mean <= hi
