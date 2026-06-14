import numpy as np
from afib_microgravity.cell_model import APCell

def test_apcell_roundtrip_shapes_and_rest():
    cell = APCell(shape=(8, 8))
    assert cell.V.shape == (8, 8)
    I = cell.reaction_step(dt=0.02)
    assert I.shape == (8, 8)
    assert np.allclose(cell.V, 0.0, atol=1e-9)  # no stimulus => stays at rest
