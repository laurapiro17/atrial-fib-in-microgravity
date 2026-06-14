import numpy as np
from afib_microgravity.crn import CRNCell, ap_biomarkers


# standard suprathreshold stimulus calibrated so V_peak is physiological (~29 mV)
def _pace_single_cell(cell, bcl=1000.0, n_beats=8, dt=0.02, stim_amp=-2300.0,
                      stim_dur=2.0):
    """Pace a 1x1 CRN cell at basic cycle length bcl (ms); return (t, V) of last beat.
    NOTE: reaction_step returns dV/dt and does NOT mutate V, so this helper updates V."""
    steps_per_beat = int(bcl / dt)
    stim_steps = int(stim_dur / dt)
    t_last, v_last = [], []
    for beat in range(n_beats):
        for s in range(steps_per_beat):
            I_stim = stim_amp if s < stim_steps else 0.0
            dVdt = cell.reaction_step(dt, I_stim=I_stim)
            cell.set_V(cell.V + dt * dVdt)
            if beat == n_beats - 1:
                v_last.append(float(cell.V.ravel()[0]))
                t_last.append(s * dt)
    return np.array(t_last), np.array(v_last)


def test_crn_single_cell_biomarkers_in_physiological_range():
    cell = CRNCell(shape=(1, 1))
    t, V = _pace_single_cell(cell)
    bm = ap_biomarkers(t, V)
    assert -86.0 <= bm["V_rest"] <= -75.0,  bm
    assert  10.0 <= bm["V_peak"] <= 35.0,   bm
    assert 250.0 <= bm["APD90"]  <= 360.0,  bm
    assert bm["dVdt_max"] >= 100.0,          bm


def test_crn_apd_shortens_at_faster_rate():
    slow = CRNCell(shape=(1, 1)); fast = CRNCell(shape=(1, 1))
    ts, Vs = _pace_single_cell(slow, bcl=1000.0)
    tf, Vf = _pace_single_cell(fast, bcl=400.0)
    assert ap_biomarkers(tf, Vf)["APD90"] < ap_biomarkers(ts, Vs)["APD90"]
