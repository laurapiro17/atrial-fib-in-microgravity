"""Numba-accelerated CRN reaction kernel.

This is an *exact* transcription of :meth:`afib_microgravity.crn.CRNCell.reaction_step`
into a single ``@njit(parallel=True)`` loop over the flattened grid. The per-cell
math (currents, Rush-Larsen gates, forward-Euler concentrations) is identical to the
pure-NumPy path, so single-cell biomarkers are unchanged; only the execution engine
differs. ``crn.py`` imports ``HAVE_NUMBA`` / ``reaction_kernel`` from here and falls
back to the NumPy implementation when Numba is unavailable.

State arrays are 1-D float64 views (``ravel()``) of the grid and are mutated in place.
The kernel returns ``dVdt`` (mV/ms) as a fresh 1-D array.
"""
from __future__ import annotations

import numpy as np

try:
    from numba import njit, prange

    HAVE_NUMBA = True
except Exception:  # pragma: no cover - exercised only when numba missing
    HAVE_NUMBA = False

    def njit(*args, **kwargs):  # type: ignore
        def wrap(f):
            return f

        if args and callable(args[0]):
            return args[0]
        return wrap

    prange = range  # type: ignore


@njit(cache=True, parallel=True, fastmath=False)
def reaction_kernel(
    V, m, h, j, oa, oi, ua, ui, xr, xs, d, f, fca,
    u_rel, v_rel, w_rel, Nai, Ki, Cai, Ca_up, Ca_rel,
    dt, I_stim,
    # scalar params
    R, T, F, Cm, Nao, Ko, Cao,
    g_Na, g_K1, g_to, g_Kur_scale, g_Kr, g_Ks, g_CaL, g_bCa, g_bNa,
    I_NaK_max, Km_Nai, Km_Ko, I_NaCa_max, K_mNa, K_mCa, K_sat, gamma,
    I_pCa_max, Km_pCa, K_rel, tau_tr, I_up_max, K_up, Ca_up_max,
    CMDN_max, TRPN_max, CSQN_max, Km_CMDN, Km_TRPN, Km_CSQN,
    Vi, Vup, Vrel, K_Q10,
):
    n = V.shape[0]
    dVdt = np.empty(n, dtype=np.float64)
    RTF = R * T / F
    sigma = (np.exp(Nao / 67.3) - 1.0) / 7.0
    FRT = F / (R * T)
    I_stim_norm = I_stim / Cm

    for k in prange(n):
        Vk = V[k]
        Naik = Nai[k]
        Kik = Ki[k]
        Caik = Cai[k]

        E_Na = RTF * np.log(Nao / Naik)
        E_K = RTF * np.log(Ko / Kik)
        E_Ca = 0.5 * RTF * np.log(Cao / Caik)

        # ---- INa ----
        denom_m = 1.0 - np.exp(-0.1 * (Vk + 47.13))
        if abs(Vk + 47.13) < 1e-10:
            alpha_m = 3.2
        else:
            alpha_m = 0.32 * (Vk + 47.13) / denom_m
        beta_m = 0.08 * np.exp(-Vk / 11.0)
        m_inf = alpha_m / (alpha_m + beta_m)
        tau_m = 1.0 / (alpha_m + beta_m)

        if Vk < -40.0:
            alpha_h = 0.135 * np.exp(-(Vk + 80.0) / 6.8)
            beta_h = 3.56 * np.exp(0.079 * Vk) + 3.1e5 * np.exp(0.35 * Vk)
            alpha_j = ((-1.2714e5 * np.exp(0.2444 * Vk)
                        - 3.474e-5 * np.exp(-0.04391 * Vk))
                       * (Vk + 37.78) / (1.0 + np.exp(0.311 * (Vk + 79.23))))
            beta_j = (0.1212 * np.exp(-0.01052 * Vk)
                      / (1.0 + np.exp(-0.1378 * (Vk + 40.14))))
        else:
            alpha_h = 0.0
            beta_h = 1.0 / (0.13 * (1.0 + np.exp(-(Vk + 10.66) / 11.1)))
            alpha_j = 0.0
            beta_j = (0.3 * np.exp(-2.535e-7 * Vk)
                      / (1.0 + np.exp(-0.1 * (Vk + 32.0))))
        h_inf = alpha_h / (alpha_h + beta_h)
        tau_h = 1.0 / (alpha_h + beta_h)
        j_inf = alpha_j / (alpha_j + beta_j)
        tau_j = 1.0 / (alpha_j + beta_j)

        mk = m[k]
        I_Na = g_Na * mk * mk * mk * h[k] * j[k] * (Vk - E_Na)

        # ---- Ito ----
        alpha_oa = 0.65 / (np.exp(-(Vk + 10.0) / 8.5)
                           + np.exp(-(Vk - 30.0) / 59.0))
        beta_oa = 0.65 / (2.5 + np.exp((Vk + 82.0) / 17.0))
        tau_oa = 1.0 / (alpha_oa + beta_oa) / K_Q10
        oa_inf = 1.0 / (1.0 + np.exp(-(Vk + 20.47) / 17.54))

        alpha_oi = 1.0 / (18.53 + np.exp((Vk + 113.7) / 10.95))
        beta_oi = 1.0 / (35.56 + np.exp(-(Vk + 1.26) / 7.44))
        tau_oi = 1.0 / (alpha_oi + beta_oi) / K_Q10
        oi_inf = 1.0 / (1.0 + np.exp((Vk + 43.1) / 5.3))

        oak = oa[k]
        I_to = g_to * oak * oak * oak * oi[k] * (Vk - E_K)

        # ---- IKur ----
        alpha_ua = 0.65 / (np.exp(-(Vk + 10.0) / 8.5)
                           + np.exp(-(Vk - 30.0) / 59.0))
        beta_ua = 0.65 / (2.5 + np.exp((Vk + 82.0) / 17.0))
        tau_ua = 1.0 / (alpha_ua + beta_ua) / K_Q10
        ua_inf = 1.0 / (1.0 + np.exp(-(Vk + 30.3) / 9.6))

        alpha_ui = 1.0 / (21.0 + np.exp(-(Vk - 185.0) / 28.0))
        beta_ui = np.exp((Vk - 158.0) / 16.0)
        tau_ui = 1.0 / (alpha_ui + beta_ui) / K_Q10
        ui_inf = 1.0 / (1.0 + np.exp((Vk - 99.45) / 27.48))

        g_Kur = g_Kur_scale * (0.005 + 0.05 / (1.0 + np.exp(-(Vk - 15.0) / 13.0)))
        uak = ua[k]
        I_Kur = g_Kur * uak * uak * uak * ui[k] * (Vk - E_K)

        # ---- IKr ----
        denom_xr_a = 1.0 - np.exp(-(Vk + 14.1) / 5.0)
        if abs(Vk + 14.1) < 1e-10:
            alpha_xr = 0.0015
        else:
            alpha_xr = 0.0003 * (Vk + 14.1) / denom_xr_a
        denom_xr_b = np.exp((Vk - 3.3328) / 5.1237) - 1.0
        if abs(Vk - 3.3328) < 1e-10:
            beta_xr = 3.7836118e-4
        else:
            beta_xr = 7.3898e-5 * (Vk - 3.3328) / denom_xr_b
        tau_xr = 1.0 / (alpha_xr + beta_xr)
        xr_inf = 1.0 / (1.0 + np.exp(-(Vk + 14.1) / 6.5))
        I_Kr = g_Kr * xr[k] * (Vk - E_K) / (1.0 + np.exp((Vk + 15.0) / 22.4))

        # ---- IKs ----
        denom_xs_a = 1.0 - np.exp(-(Vk - 19.9) / 17.0)
        if abs(Vk - 19.9) < 1e-10:
            alpha_xs = 0.00068
        else:
            alpha_xs = 4e-5 * (Vk - 19.9) / denom_xs_a
        denom_xs_b = np.exp((Vk - 19.9) / 9.0) - 1.0
        if abs(Vk - 19.9) < 1e-10:
            beta_xs = 0.000315
        else:
            beta_xs = 3.5e-5 * (Vk - 19.9) / denom_xs_b
        tau_xs = 0.5 / (alpha_xs + beta_xs)
        xs_inf = 1.0 / np.sqrt(1.0 + np.exp(-(Vk - 19.9) / 12.7))
        xsk = xs[k]
        I_Ks = g_Ks * xsk * xsk * (Vk - E_K)

        # ---- ICaL ----
        d_inf = 1.0 / (1.0 + np.exp(-(Vk + 10.0) / 8.0))
        if abs(Vk + 10.0) < 1e-10:
            tau_d = 4.579 / (1.0 + np.exp(-(Vk + 10.0) / 6.24))
        else:
            tau_d = ((1.0 - np.exp(-(Vk + 10.0) / 6.24))
                     / (0.035 * (Vk + 10.0)
                        * (1.0 + np.exp(-(Vk + 10.0) / 6.24))))
        f_inf = 1.0 / (1.0 + np.exp((Vk + 28.0) / 6.9))
        tau_f = 9.0 / (0.0197 * np.exp(-(0.0337 * 0.0337)
                                       * (Vk + 10.0) * (Vk + 10.0)) + 0.02)
        fca_inf = 1.0 / (1.0 + Caik / 0.00035)
        tau_fca = 2.0
        I_CaL = g_CaL * d[k] * f[k] * fca[k] * (Vk - 65.0)

        # ---- INaK ----
        f_NaK = 1.0 / (1.0 + 0.1245 * np.exp(-0.1 * FRT * Vk)
                       + 0.0365 * sigma * np.exp(-FRT * Vk))
        I_NaK = (I_NaK_max * f_NaK
                 * (1.0 / (1.0 + (Km_Nai / Naik) ** 1.5))
                 * Ko / (Ko + Km_Ko))

        # ---- INaCa ----
        exp_g = np.exp(gamma * FRT * Vk)
        exp_g1 = np.exp((gamma - 1.0) * FRT * Vk)
        I_NaCa = (I_NaCa_max
                  * (exp_g * Naik * Naik * Naik * Cao
                     - exp_g1 * Nao * Nao * Nao * Caik)
                  / ((K_mNa ** 3 + Nao ** 3) * (K_mCa + Cao)
                     * (1.0 + K_sat * exp_g1)))

        # ---- background ----
        I_bCa = g_bCa * (Vk - E_Ca)
        I_bNa = g_bNa * (Vk - E_Na)

        # ---- IpCa ----
        I_pCa = I_pCa_max * Caik / (Km_pCa + Caik)

        # ---- IK1 ----
        I_K1 = g_K1 * (Vk - E_K) / (1.0 + np.exp(0.07 * (Vk + 80.0)))

        # ---- SR fluxes ----
        I_rel = K_rel * u_rel[k] * u_rel[k] * v_rel[k] * w_rel[k] * (Ca_rel[k] - Caik)
        I_tr = (Ca_up[k] - Ca_rel[k]) / tau_tr
        I_up = I_up_max / (1.0 + K_up / Caik)
        I_up_leak = I_up_max * Ca_up[k] / Ca_up_max

        I_ion = (I_Na + I_K1 + I_to + I_Kur + I_Kr + I_Ks
                 + I_CaL + I_NaK + I_NaCa + I_bNa + I_bCa + I_pCa)
        dVdt[k] = -(I_ion + I_stim_norm)

        # ---- Rush-Larsen gates ----
        m[k] = m_inf - (m_inf - mk) * np.exp(-dt / tau_m)
        h[k] = h_inf - (h_inf - h[k]) * np.exp(-dt / tau_h)
        j[k] = j_inf - (j_inf - j[k]) * np.exp(-dt / tau_j)
        oa[k] = oa_inf - (oa_inf - oak) * np.exp(-dt / tau_oa)
        oi[k] = oi_inf - (oi_inf - oi[k]) * np.exp(-dt / tau_oi)
        ua[k] = ua_inf - (ua_inf - uak) * np.exp(-dt / tau_ua)
        ui[k] = ui_inf - (ui_inf - ui[k]) * np.exp(-dt / tau_ui)
        xr[k] = xr_inf - (xr_inf - xr[k]) * np.exp(-dt / tau_xr)
        xs[k] = xs_inf - (xs_inf - xsk) * np.exp(-dt / tau_xs)
        d[k] = d_inf - (d_inf - d[k]) * np.exp(-dt / tau_d)
        f[k] = f_inf - (f_inf - f[k]) * np.exp(-dt / tau_f)
        fca[k] = fca_inf - (fca_inf - fca[k]) * np.exp(-dt / tau_fca)

        # ---- SR release gates ----
        Fn = (1e-12 * Vrel * I_rel
              - (5e-13 / F) * (0.5 * I_CaL * Cm - 0.2 * I_NaCa * Cm))
        u_inf = 1.0 / (1.0 + np.exp(-(Fn - 3.4175e-13) / 13.67e-16))
        tau_u = 8.0
        v_inf = 1.0 - 1.0 / (1.0 + np.exp(-(Fn - 6.835e-14) / 13.67e-16))
        tau_v = 1.91 + 2.09 / (1.0 + np.exp(-(Fn - 3.4175e-13) / 13.67e-16))
        w_inf = 1.0 - 1.0 / (1.0 + np.exp(-(Vk - 40.0) / 17.0))
        if abs(Vk - 7.9) < 1e-10:
            tau_w = 6.0 * 0.2 / 1.3
        else:
            tau_w = (6.0 * (1.0 - np.exp(-(Vk - 7.9) / 5.0))
                     / ((1.0 + 0.3 * np.exp(-(Vk - 7.9) / 5.0)) * (Vk - 7.9)))
        u_rel[k] = u_inf - (u_inf - u_rel[k]) * np.exp(-dt / tau_u)
        v_rel[k] = v_inf - (v_inf - v_rel[k]) * np.exp(-dt / tau_v)
        w_rel[k] = w_inf - (w_inf - w_rel[k]) * np.exp(-dt / tau_w)

        # ---- concentrations (forward Euler) ----
        dNai = (-3.0 * I_NaK - 3.0 * I_NaCa - I_bNa - I_Na) * Cm / (F * Vi)
        dKi = (2.0 * I_NaK - I_K1 - I_to - I_Kur - I_Kr - I_Ks) * Cm / (F * Vi)
        B1 = ((2.0 * I_NaCa - I_pCa - I_CaL - I_bCa) * Cm / (2.0 * F * Vi)
              + (Vup * (I_up_leak - I_up) + I_rel * Vrel) / Vi)
        B2 = 1.0 + (TRPN_max * Km_TRPN / (Caik + Km_TRPN) ** 2
                    + CMDN_max * Km_CMDN / (Caik + Km_CMDN) ** 2)
        dCai = B1 / B2
        dCa_up = I_up - I_up_leak - I_tr * Vrel / Vup
        dCa_rel = ((I_tr - I_rel)
                   * (1.0 / (1.0 + CSQN_max * Km_CSQN
                             / (Ca_rel[k] + Km_CSQN) ** 2)))

        Nai[k] = Naik + dt * dNai
        Ki[k] = Kik + dt * dKi
        new_Cai = Caik + dt * dCai
        Cai[k] = new_Cai if new_Cai > 1e-7 else 1e-7
        Ca_up[k] = Ca_up[k] + dt * dCa_up
        Ca_rel[k] = Ca_rel[k] + dt * dCa_rel

    return dVdt
