"""APD90 for ground vs microgravity (severity sweep) + wavelength + domain size."""
import os, sys, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
from afib_microgravity.crn import CRNCell, ap_biomarkers
from afib_microgravity.remodeling import baseline_crn_params, microgravity_crn_params

DT = 0.02
D_LONG = 0.15  # chosen: planar CV ~58 cm/s


def apd90(params, bcl=1000.0, n_beats=12):
    cell = CRNCell(shape=(1, 1), params=params, use_numba=True)
    spb = int(bcl / DT); ss = int(2.0 / DT)
    tl, vl = [], []
    for beat in range(n_beats):
        for s in range(spb):
            I = -2300.0 if s < ss else 0.0
            dV = cell.reaction_step(DT, I_stim=I)
            cell.set_V(cell.V + DT * dV)
            if beat == n_beats - 1:
                vl.append(float(cell.V.ravel()[0])); tl.append(s * DT)
    return ap_biomarkers(np.array(tl), np.array(vl))


if __name__ == "__main__":
    # CV ~ const for fixed d_long & unchanged g_Na; reuse measured 58 cm/s.
    CV = 58.0
    print(f"d_long={D_LONG}  planar CV~{CV} cm/s (g_Na unchanged across conditions)\n", flush=True)

    print("=== Ground (baseline) ===", flush=True)
    bg = apd90(baseline_crn_params())
    print(f"  APD90={bg['APD90']:.1f} ms  Vrest={bg['V_rest']:.1f}  Vpeak={bg['V_peak']:.1f}  dVdt={bg['dVdt_max']:.0f}", flush=True)
    wl_g = CV * (bg['APD90'] / 1000.0)  # cm/s * s = cm
    print(f"  wavelength = CV*APD = {wl_g:.2f} cm  -> rotor needs ~{wl_g:.1f}cm = {wl_g/0.025:.0f} cells/side\n", flush=True)

    for sev in (0.5, 0.75, 1.0):
        print(f"=== Microgravity severity={sev} ===", flush=True)
        bm = apd90(microgravity_crn_params(severity=sev))
        print(f"  APD90={bm['APD90']:.1f} ms  Vrest={bm['V_rest']:.1f}  Vpeak={bm['V_peak']:.1f}  dVdt={bm['dVdt_max']:.0f}", flush=True)
        wl = CV * (bm['APD90'] / 1000.0)
        print(f"  wavelength = {wl:.2f} cm  -> rotor needs ~{wl:.1f}cm = {wl/0.025:.0f} cells/side\n", flush=True)
