import numpy as np
from afib_microgravity.diffusion import AnisotropicDiffusion

def test_isotropic_limit_matches_5point_laplacian():
    rng = np.random.default_rng(0)
    u = rng.standard_normal((16, 16))
    op = AnisotropicDiffusion(shape=(16,16), d_long=1.0, d_trans=1.0,
                              theta=np.zeros((16,16)), dx=1.0)
    lap = (np.roll(u,1,0)+np.roll(u,-1,0)+np.roll(u,1,1)+np.roll(u,-1,1)-4*u)
    assert np.allclose(op.apply(u)[1:-1,1:-1], lap[1:-1,1:-1], atol=1e-9)

def test_no_flux_conserves_total_for_uniform_field():
    u = np.ones((10, 10)) * 0.7
    op = AnisotropicDiffusion(shape=(10,10), d_long=2.0, d_trans=0.5,
                              theta=np.full((10,10), 0.3), dx=0.25)
    assert np.allclose(op.apply(u), 0.0, atol=1e-9)
