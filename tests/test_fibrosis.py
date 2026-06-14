import numpy as np
from afib_microgravity.fibrosis import correlated_fibrosis

def test_density_is_approximately_hit():
    f = correlated_fibrosis((128,128), density=0.3, corr_len=4.0,
                            d_healthy=1.0, d_scar=0.05, seed=1)
    frac_scar = np.mean(f < 0.5)
    assert abs(frac_scar - 0.3) < 0.05

def test_zero_density_is_homogeneous():
    f = correlated_fibrosis((32,32), density=0.0, corr_len=4.0, seed=1)
    assert np.allclose(f, 1.0)

def test_larger_corr_len_makes_larger_patches():
    small = correlated_fibrosis((128,128), 0.3, corr_len=1.0, seed=2)
    large = correlated_fibrosis((128,128), 0.3, corr_len=8.0, seed=2)
    def autocorr(a, lag):
        b = (a < 0.5).astype(float); b -= b.mean()
        return np.mean(b * np.roll(b, lag, axis=0))
    assert autocorr(large, 3) > autocorr(small, 3)
