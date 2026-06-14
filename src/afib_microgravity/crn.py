"""Courtemanche-Ramirez-Nattel (1998) human atrial single-cell ionic model.

Vectorised over an arbitrary grid: every state variable is a 2D ndarray of the
cell's ``shape`` and every current is an array operation, so the same kinetics
drive a 1x1 single cell or a (200, 200) monodomain sheet without Python loops.

Gates are advanced with Rush-Larsen; V and concentrations with forward Euler.
Equations and constants transcribed from the CellML listing of the model in the
Physiome Model Repository (courtemanche_ramirez_nattel_1998.cellml), which
reproduces Courtemanche, Ramirez & Nattel, Am J Physiol 1998;275:H301-H321.

Currents are in pA/pF (Cm-normalised), so the membrane derivative returned by
``reaction_step`` is dV/dt = -(I_ion_total + I_stim) in mV/ms.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .cell_model import CellModel
from .crn_numba import HAVE_NUMBA, reaction_kernel

# Names of the 21 mutable state arrays, in the order the Numba kernel expects.
_STATE_NAMES = (
    "m", "h", "j", "oa", "oi", "ua", "ui", "xr", "xs", "d", "f", "fca",
    "u_rel", "v_rel", "w_rel", "Nai", "Ki", "Cai", "Ca_up", "Ca_rel",
)


@dataclass
class CRNParams:
    # Physical constants
    R: float = 8.3143          # J / (mol K)
    T: float = 310.0           # K
    F: float = 96.4867         # C / mmol
    Cm: float = 100.0          # pF

    # External concentrations (mM)
    Nao: float = 140.0
    Ko: float = 5.4
    Cao: float = 1.8

    # Maximal conductances / scalings
    g_Na: float = 7.8          # nS/pF
    g_K1: float = 0.09
    g_to: float = 0.1652
    g_Kur_scale: float = 1.0   # multiplies the voltage-dependent g_Kur
    g_Kr: float = 0.029411765
    g_Ks: float = 0.12941176
    g_CaL: float = 0.12375
    g_bCa: float = 0.001131
    g_bNa: float = 0.0006744375

    # Pumps and exchanger
    I_NaK_max: float = 0.59933874   # pA/pF
    Km_Nai: float = 10.0            # mM
    Km_Ko: float = 1.5              # mM
    I_NaCa_max: float = 1600.0      # pA/pF
    K_mNa: float = 87.5             # mM
    K_mCa: float = 1.38             # mM
    K_sat: float = 0.1
    gamma: float = 0.35
    I_pCa_max: float = 0.275        # pA/pF
    Km_pCa: float = 0.0005          # mM (Ca half-sat for IpCa, = 0.5 uM)

    # SR Ca handling
    K_rel: float = 30.0             # 1/ms
    tau_tr: float = 180.0           # ms
    I_up_max: float = 0.005         # mM/ms
    K_up: float = 0.00092           # mM
    Ca_up_max: float = 15.0         # mM

    # Buffering (mM)
    CMDN_max: float = 0.05
    TRPN_max: float = 0.07
    CSQN_max: float = 10.0
    Km_CMDN: float = 0.00238
    Km_TRPN: float = 0.0005
    Km_CSQN: float = 0.8

    # Cell geometry (uL -> for concentration ODEs)
    Vcell: float = 20100.0          # um^3
    Vi: float = 13668.0             # um^3 (0.68 * Vcell)
    Vup: float = 1109.52            # um^3 (0.0552 * Vcell)
    Vrel: float = 96.48             # um^3 (0.0048 * Vcell)

    K_Q10: float = 3.0


# Published CRN 1998 resting initial conditions (Table / appendix values).
_IC = dict(
    V=-81.18, m=2.908e-3, h=9.649e-1, j=9.775e-1,
    oa=3.043e-2, oi=9.992e-1, ua=4.966e-3, ui=9.986e-1,
    xr=3.296e-5, xs=1.869e-2, d=1.367e-4, f=9.996e-1, fca=7.755e-1,
    u_rel=0.0, v_rel=1.0, w_rel=9.992e-1,
    Nai=11.17, Ki=139.0, Cai=1.013e-4, Ca_up=1.488, Ca_rel=1.488,
)


def _full(shape, value):
    return np.full(shape, value, dtype=np.float64)


class CRNCell(CellModel):
    """CRN 1998 human-atrial kinetics behind the CellModel interface."""

    def __init__(self, shape=(200, 200), params: CRNParams | None = None,
                 use_numba: bool | None = None):
        self.p = params or CRNParams()
        self.shape = shape
        # use_numba=None -> auto (on when numba importable); explicit overrides.
        self.use_numba = HAVE_NUMBA if use_numba is None else bool(use_numba)
        for name, val in _IC.items():
            if name == "V":
                self._Vm = _full(shape, val)
            else:
                setattr(self, name, _full(shape, val))

    # ---- CellModel interface ------------------------------------------------
    @property
    def V(self) -> np.ndarray:
        return self._Vm

    def set_V(self, V: np.ndarray) -> None:
        self._Vm = V

    def stimulate(self, mask: np.ndarray, value: float) -> None:
        self._Vm[mask] = value

    # ---- kinetics -----------------------------------------------------------
    def reaction_step(self, dt: float, I_stim: float = 0.0) -> np.ndarray:
        if self.use_numba:
            return self._reaction_step_numba(dt, I_stim)
        return self._reaction_step_numpy(dt, I_stim)

    def _reaction_step_numba(self, dt: float, I_stim: float = 0.0) -> np.ndarray:
        """Dispatch the per-cell kinetics to the compiled Numba kernel.

        State arrays are mutated in place through flat (ravel) views; only dV/dt
        is returned (reshaped to the grid), matching the NumPy path's contract.
        """
        p = self.p
        flat = [self._Vm.reshape(-1)] + [
            getattr(self, n).reshape(-1) for n in _STATE_NAMES
        ]
        dVdt = reaction_kernel(
            *flat, float(dt), float(I_stim),
            p.R, p.T, p.F, p.Cm, p.Nao, p.Ko, p.Cao,
            p.g_Na, p.g_K1, p.g_to, p.g_Kur_scale, p.g_Kr, p.g_Ks, p.g_CaL,
            p.g_bCa, p.g_bNa,
            p.I_NaK_max, p.Km_Nai, p.Km_Ko, p.I_NaCa_max, p.K_mNa, p.K_mCa,
            p.K_sat, p.gamma, p.I_pCa_max, p.Km_pCa, p.K_rel, p.tau_tr,
            p.I_up_max, p.K_up, p.Ca_up_max,
            p.CMDN_max, p.TRPN_max, p.CSQN_max, p.Km_CMDN, p.Km_TRPN, p.Km_CSQN,
            p.Vi, p.Vup, p.Vrel, p.K_Q10,
        )
        return dVdt.reshape(self.shape)

    def _reaction_step_numpy(self, dt: float, I_stim: float = 0.0) -> np.ndarray:
        p = self.p
        V = self._Vm
        RTF = p.R * p.T / p.F

        # Nernst potentials
        E_Na = RTF * np.log(p.Nao / self.Nai)
        E_K = RTF * np.log(p.Ko / self.Ki)
        E_Ca = 0.5 * RTF * np.log(p.Cao / self.Cai)

        # ---- INa ----
        denom_m = 1.0 - np.exp(-0.1 * (V + 47.13))
        alpha_m = np.where(np.abs(V + 47.13) < 1e-10,
                           3.2,
                           0.32 * (V + 47.13) / denom_m)
        beta_m = 0.08 * np.exp(-V / 11.0)
        m_inf = alpha_m / (alpha_m + beta_m)
        tau_m = 1.0 / (alpha_m + beta_m)

        lo = V < -40.0
        alpha_h = np.where(lo, 0.135 * np.exp(-(V + 80.0) / 6.8), 0.0)
        beta_h = np.where(
            lo,
            3.56 * np.exp(0.079 * V) + 3.1e5 * np.exp(0.35 * V),
            1.0 / (0.13 * (1.0 + np.exp(-(V + 10.66) / 11.1))),
        )
        h_inf = alpha_h / (alpha_h + beta_h)
        tau_h = 1.0 / (alpha_h + beta_h)

        alpha_j = np.where(
            lo,
            (-1.2714e5 * np.exp(0.2444 * V) - 3.474e-5 * np.exp(-0.04391 * V))
            * (V + 37.78) / (1.0 + np.exp(0.311 * (V + 79.23))),
            0.0,
        )
        beta_j = np.where(
            lo,
            0.1212 * np.exp(-0.01052 * V) / (1.0 + np.exp(-0.1378 * (V + 40.14))),
            0.3 * np.exp(-2.535e-7 * V) / (1.0 + np.exp(-0.1 * (V + 32.0))),
        )
        j_inf = alpha_j / (alpha_j + beta_j)
        tau_j = 1.0 / (alpha_j + beta_j)

        I_Na = p.g_Na * self.m**3 * self.h * self.j * (V - E_Na)

        # ---- Ito ----
        K_Q10 = p.K_Q10
        alpha_oa = 0.65 / (np.exp(-(V + 10.0) / 8.5) + np.exp(-(V - 30.0) / 59.0))
        beta_oa = 0.65 / (2.5 + np.exp((V + 82.0) / 17.0))
        tau_oa = 1.0 / (alpha_oa + beta_oa) / K_Q10
        oa_inf = 1.0 / (1.0 + np.exp(-(V + 20.47) / 17.54))

        alpha_oi = 1.0 / (18.53 + np.exp((V + 113.7) / 10.95))
        beta_oi = 1.0 / (35.56 + np.exp(-(V + 1.26) / 7.44))
        tau_oi = 1.0 / (alpha_oi + beta_oi) / K_Q10
        oi_inf = 1.0 / (1.0 + np.exp((V + 43.1) / 5.3))

        I_to = p.g_to * self.oa**3 * self.oi * (V - E_K)

        # ---- IKur ----
        alpha_ua = 0.65 / (np.exp(-(V + 10.0) / 8.5) + np.exp(-(V - 30.0) / 59.0))
        beta_ua = 0.65 / (2.5 + np.exp((V + 82.0) / 17.0))
        tau_ua = 1.0 / (alpha_ua + beta_ua) / K_Q10
        ua_inf = 1.0 / (1.0 + np.exp(-(V + 30.3) / 9.6))

        alpha_ui = 1.0 / (21.0 + np.exp(-(V - 185.0) / 28.0))
        beta_ui = np.exp((V - 158.0) / 16.0)
        tau_ui = 1.0 / (alpha_ui + beta_ui) / K_Q10
        ui_inf = 1.0 / (1.0 + np.exp((V - 99.45) / 27.48))

        g_Kur = p.g_Kur_scale * (0.005 + 0.05 / (1.0 + np.exp(-(V - 15.0) / 13.0)))
        I_Kur = g_Kur * self.ua**3 * self.ui * (V - E_K)

        # ---- IKr ----
        denom_xr_a = 1.0 - np.exp(-(V + 14.1) / 5.0)
        alpha_xr = np.where(np.abs(V + 14.1) < 1e-10,
                            0.0015,
                            0.0003 * (V + 14.1) / denom_xr_a)
        denom_xr_b = np.exp((V - 3.3328) / 5.1237) - 1.0
        beta_xr = np.where(np.abs(V - 3.3328) < 1e-10,
                           3.7836118e-4,
                           7.3898e-5 * (V - 3.3328) / denom_xr_b)
        tau_xr = 1.0 / (alpha_xr + beta_xr)
        xr_inf = 1.0 / (1.0 + np.exp(-(V + 14.1) / 6.5))
        I_Kr = p.g_Kr * self.xr * (V - E_K) / (1.0 + np.exp((V + 15.0) / 22.4))

        # ---- IKs ----
        denom_xs_a = 1.0 - np.exp(-(V - 19.9) / 17.0)
        alpha_xs = np.where(np.abs(V - 19.9) < 1e-10,
                            0.00068,
                            4e-5 * (V - 19.9) / denom_xs_a)
        denom_xs_b = np.exp((V - 19.9) / 9.0) - 1.0
        beta_xs = np.where(np.abs(V - 19.9) < 1e-10,
                           0.000315,
                           3.5e-5 * (V - 19.9) / denom_xs_b)
        tau_xs = 0.5 / (alpha_xs + beta_xs)
        xs_inf = 1.0 / np.sqrt(1.0 + np.exp(-(V - 19.9) / 12.7))
        I_Ks = p.g_Ks * self.xs**2 * (V - E_K)

        # ---- ICaL ----
        d_inf = 1.0 / (1.0 + np.exp(-(V + 10.0) / 8.0))
        denom_d = 1.0 - np.exp(-(V + 10.0) / 6.24)
        tau_d = np.where(
            np.abs(V + 10.0) < 1e-10,
            4.579 / (1.0 + np.exp(-(V + 10.0) / 6.24)),
            (1.0 - np.exp(-(V + 10.0) / 6.24))
            / (0.035 * (V + 10.0) * (1.0 + np.exp(-(V + 10.0) / 6.24))),
        )
        _ = denom_d
        f_inf = 1.0 / (1.0 + np.exp((V + 28.0) / 6.9))
        tau_f = 9.0 / (0.0197 * np.exp(-(0.0337**2) * (V + 10.0)**2) + 0.02)
        fca_inf = 1.0 / (1.0 + self.Cai / 0.00035)
        tau_fca = 2.0  # ms
        I_CaL = p.g_CaL * self.d * self.f * self.fca * (V - 65.0)

        # ---- INaK ----
        sigma = (np.exp(p.Nao / 67.3) - 1.0) / 7.0
        f_NaK = 1.0 / (1.0 + 0.1245 * np.exp(-0.1 * p.F * V / (p.R * p.T))
                       + 0.0365 * sigma * np.exp(-p.F * V / (p.R * p.T)))
        I_NaK = (p.I_NaK_max * f_NaK
                 * (1.0 / (1.0 + (p.Km_Nai / self.Nai) ** 1.5))
                 * p.Ko / (p.Ko + p.Km_Ko))

        # ---- INaCa ----
        exp_g = np.exp(p.gamma * p.F * V / (p.R * p.T))
        exp_g1 = np.exp((p.gamma - 1.0) * p.F * V / (p.R * p.T))
        I_NaCa = (
            p.I_NaCa_max
            * (exp_g * self.Nai**3 * p.Cao - exp_g1 * p.Nao**3 * self.Cai)
            / ((p.K_mNa**3 + p.Nao**3) * (p.K_mCa + p.Cao)
               * (1.0 + p.K_sat * exp_g1))
        )

        # ---- background currents ----
        I_bCa = p.g_bCa * (V - E_Ca)
        I_bNa = p.g_bNa * (V - E_Na)

        # ---- IpCa ----
        I_pCa = p.I_pCa_max * self.Cai / (p.Km_pCa + self.Cai)

        # ---- SR fluxes ----
        I_rel = p.K_rel * self.u_rel**2 * self.v_rel * self.w_rel * (self.Ca_rel - self.Cai)
        I_tr = (self.Ca_up - self.Ca_rel) / p.tau_tr
        I_up = p.I_up_max / (1.0 + p.K_up / self.Cai)
        I_up_leak = p.I_up_max * self.Ca_up / p.Ca_up_max

        # ---- total ionic current (pA/pF) ----
        I_ion = (I_Na + I_K1(V, E_K, p) + I_to + I_Kur + I_Kr + I_Ks
                 + I_CaL + I_NaK + I_NaCa + I_bNa + I_bCa + I_pCa)

        # I_stim is supplied as a total stimulus (pA); normalise by Cm to pA/pF.
        dVdt = -(I_ion + I_stim / p.Cm)

        # ---- Rush-Larsen gate updates (in place) ----
        self.m = m_inf - (m_inf - self.m) * np.exp(-dt / tau_m)
        self.h = h_inf - (h_inf - self.h) * np.exp(-dt / tau_h)
        self.j = j_inf - (j_inf - self.j) * np.exp(-dt / tau_j)
        self.oa = oa_inf - (oa_inf - self.oa) * np.exp(-dt / tau_oa)
        self.oi = oi_inf - (oi_inf - self.oi) * np.exp(-dt / tau_oi)
        self.ua = ua_inf - (ua_inf - self.ua) * np.exp(-dt / tau_ua)
        self.ui = ui_inf - (ui_inf - self.ui) * np.exp(-dt / tau_ui)
        self.xr = xr_inf - (xr_inf - self.xr) * np.exp(-dt / tau_xr)
        self.xs = xs_inf - (xs_inf - self.xs) * np.exp(-dt / tau_xs)
        self.d = d_inf - (d_inf - self.d) * np.exp(-dt / tau_d)
        self.f = f_inf - (f_inf - self.f) * np.exp(-dt / tau_f)
        self.fca = fca_inf - (fca_inf - self.fca) * np.exp(-dt / tau_fca)

        # ---- SR release gates (Rush-Larsen) ----
        # Fn drives the activation of release
        Fn = 1e-12 * p.Vrel * I_rel - (5e-13 / p.F) * (0.5 * I_CaL * p.Cm
                                                       - 0.2 * I_NaCa * p.Cm)
        u_inf = 1.0 / (1.0 + np.exp(-(Fn - 3.4175e-13) / 13.67e-16))
        tau_u = 8.0
        v_inf = 1.0 - 1.0 / (1.0 + np.exp(-(Fn - 6.835e-14) / 13.67e-16))
        tau_v = 1.91 + 2.09 / (1.0 + np.exp(-(Fn - 3.4175e-13) / 13.67e-16))
        w_inf = 1.0 - 1.0 / (1.0 + np.exp(-(V - 40.0) / 17.0))
        tau_w = np.where(
            np.abs(V - 7.9) < 1e-10,
            6.0 * 0.2 / 1.3,
            6.0 * (1.0 - np.exp(-(V - 7.9) / 5.0))
            / ((1.0 + 0.3 * np.exp(-(V - 7.9) / 5.0)) * (V - 7.9)),
        )
        self.u_rel = u_inf - (u_inf - self.u_rel) * np.exp(-dt / tau_u)
        self.v_rel = v_inf - (v_inf - self.v_rel) * np.exp(-dt / tau_v)
        self.w_rel = w_inf - (w_inf - self.w_rel) * np.exp(-dt / tau_w)

        # ---- concentration updates (forward Euler) ----
        # unit factor: 1/(z F Vi) converts pA -> mM/ms; currents are pA/pF so
        # multiply by Cm (pF) to get pA. 2.0067e-4 etc. follow from geometry.
        F = p.F
        Vi = p.Vi
        Vup = p.Vup
        Vrel = p.Vrel
        Cm = p.Cm

        dNai = (-3.0 * I_NaK - 3.0 * I_NaCa - I_bNa - I_Na) * Cm / (F * Vi)
        dKi = (2.0 * I_NaK - I_K1(V, E_K, p) - I_to - I_Kur - I_Kr - I_Ks
               - I_bK(V, E_K)) * Cm / (F * Vi)

        B1 = ((2.0 * I_NaCa - I_pCa - I_CaL - I_bCa) * Cm / (2.0 * F * Vi)
              + (Vup * (I_up_leak - I_up) + I_rel * Vrel) / Vi)
        B2 = 1.0 + (p.TRPN_max * p.Km_TRPN / (self.Cai + p.Km_TRPN) ** 2
                    + p.CMDN_max * p.Km_CMDN / (self.Cai + p.Km_CMDN) ** 2)
        dCai = B1 / B2

        dCa_up = I_up - I_up_leak - I_tr * Vrel / Vup
        dCa_rel = (I_tr - I_rel) * (1.0 / (1.0 + p.CSQN_max * p.Km_CSQN
                                           / (self.Ca_rel + p.Km_CSQN) ** 2))

        self.Nai = self.Nai + dt * dNai
        self.Ki = self.Ki + dt * dKi
        self.Cai = np.maximum(self.Cai + dt * dCai, 1e-7)
        self.Ca_up = self.Ca_up + dt * dCa_up
        self.Ca_rel = self.Ca_rel + dt * dCa_rel

        return dVdt


def I_K1(V, E_K, p):
    return p.g_K1 * (V - E_K) / (1.0 + np.exp(0.07 * (V + 80.0)))


def I_bK(V, E_K):
    # CRN has no explicit background K current; placeholder kept at zero.
    return 0.0


def ap_biomarkers(t, V) -> dict:
    """Biomarkers of a single action potential trace (t in ms, V in mV)."""
    t = np.asarray(t, dtype=float)
    V = np.asarray(V, dtype=float)
    V_rest = float(np.min(V))
    V_peak = float(np.max(V))
    dVdt = np.diff(V) / np.diff(t)
    dVdt_max = float(np.max(dVdt))
    i_up = int(np.argmax(dVdt))  # index of fastest upstroke
    t_up = t[i_up]

    thr = V_rest + 0.1 * (V_peak - V_rest)  # 90% repolarisation level
    apd90 = float("nan")
    for k in range(i_up + 1, len(V)):
        if V[k] <= thr:
            # linear interpolation of the downward crossing
            v0, v1 = V[k - 1], V[k]
            t0, t1 = t[k - 1], t[k]
            t_cross = t0 + (thr - v0) * (t1 - t0) / (v1 - v0)
            apd90 = float(t_cross - t_up)
            break
    return {
        "V_rest": V_rest,
        "V_peak": V_peak,
        "dVdt_max": dVdt_max,
        "APD90": apd90,
    }
