"""APD and CV restitution (S1-S2) with the CRN human-atrial ionic model.

The restitution curve is the workhorse of arrhythmia electrophysiology: it
relates the duration of an action potential (APD) -- and the speed of the
propagating wavefront (CV) -- to the diastolic interval (DI) that preceded it.
A steep APD restitution slope promotes alternans and wavebreak, the substrate
of fibrillation, so quantifying it under a remodelled (e.g. microgravity)
parameter set is the natural single-cell read-out.

Protocol: pace ``n_s1`` conditioning S1 beats at a fixed basic cycle length to
reach a near-steady state, then deliver one premature S2 stimulus a chosen DI
after the last S1 has fully repolarised, and measure the S2 action potential.
"""
from __future__ import annotations

import numpy as np

from .crn import CRNCell, CRNParams, ap_biomarkers
from .diffusion import AnisotropicDiffusion
from .model import MonodomainSheet

_REPOL_MV = -70.0  # membrane voltage taken as "repolarised" (end of an AP)


def _scalar(V):
    """Extract the scalar voltage of a single-cell (1, 1) state."""
    return float(np.asarray(V).flat[0])


def _advance(cell, dt, n_steps, i_stim=0.0):
    """Integrate a (single-cell) CRNCell for ``n_steps`` with constant I_stim."""
    for _ in range(int(n_steps)):
        dVdt = cell.reaction_step(dt, I_stim=i_stim)
        cell.set_V(cell.V + dt * dVdt)


def apd_restitution(diastolic_intervals, dt=0.05, n_s1=4, s1_bcl=500.0,
                    stim_amp=-2300.0, stim_dur=2.0):
    """S1-S2 APD90 restitution for a single CRN atrial cell.

    For each DI in ``diastolic_intervals`` a fresh cell is paced with ``n_s1``
    S1 beats at ``s1_bcl``; an S2 is then delivered ``DI`` ms after the last S1
    has repolarised, and the S2 beat's APD90 is measured.

    Returns ``(list(diastolic_intervals), [apd90_per_di])``.
    """
    dis = list(diastolic_intervals)
    apds = []
    for di in dis:
        cell = CRNCell(shape=(1, 1), params=CRNParams())

        # --- S1 conditioning train at fixed BCL ---
        for _ in range(n_s1):
            _advance(cell, dt, round(stim_dur / dt), i_stim=stim_amp)
            _advance(cell, dt, round((s1_bcl - stim_dur) / dt), i_stim=0.0)

        # --- run out the last S1 to full repolarisation, then wait DI ms ---
        n_max = round(600.0 / dt)
        for _ in range(n_max):
            dVdt = cell.reaction_step(dt, I_stim=0.0)
            cell.set_V(cell.V + dt * dVdt)
            if _scalar(cell.V) < _REPOL_MV:
                break
        _advance(cell, dt, round(di / dt), i_stim=0.0)

        # --- S2 stimulus + record the premature beat ---
        ts = [0.0]
        vs = [_scalar(cell.V)]
        t = 0.0
        n_stim = round(stim_dur / dt)
        for k in range(n_stim):
            dVdt = cell.reaction_step(dt, I_stim=stim_amp)
            cell.set_V(cell.V + dt * dVdt)
            t += dt
            ts.append(t)
            vs.append(_scalar(cell.V))
        # record a long window to capture full S2 repolarisation
        n_record = round(500.0 / dt)
        for k in range(n_record):
            dVdt = cell.reaction_step(dt, I_stim=0.0)
            cell.set_V(cell.V + dt * dVdt)
            t += dt
            ts.append(t)
            vs.append(_scalar(cell.V))

        bm = ap_biomarkers(ts, vs)
        apds.append(bm["APD90"])
    return dis, apds


def _upstroke_time(times, volts, thr=-40.0):
    """First time the trace crosses ``thr`` upward, linearly interpolated."""
    volts = np.asarray(volts)
    times = np.asarray(times)
    for k in range(1, len(volts)):
        if volts[k - 1] < thr <= volts[k]:
            v0, v1 = volts[k - 1], volts[k]
            t0, t1 = times[k - 1], times[k]
            return t0 + (thr - v0) * (t1 - t0) / (v1 - v0)
    return float("nan")


def cv_restitution(diastolic_intervals, dt=0.05, n=60, dx=0.25, d_long=0.06,
                   d_trans=0.06, n_s1=2, s1_bcl=400.0, stim_amp=-2300.0,
                   stim_dur=2.0):
    """S1-S2 conduction-velocity restitution on a 1D CRN strand.

    A ``(1, n)`` strand is paced from its left end (S1 train then a premature
    S2 at each DI). The S2 wavefront (V crossing -40 mV upward) is timed at two
    interior probes and CV is reported in cm/s.

    Returns ``(list(diastolic_intervals), [cv_cm_per_s_per_di])``.
    """
    dis = list(diastolic_intervals)
    p1, p2 = n // 4, 3 * n // 4
    dist_cm = (p2 - p1) * dx * 0.1  # dx in mm -> cm
    cvs = []

    left = np.zeros((1, n), dtype=bool)
    left[0, 0:3] = True

    for di in dis:
        cell = CRNCell(shape=(1, n), params=CRNParams())
        diff = AnisotropicDiffusion((1, n), d_long, d_trans, theta=0.0, dx=dx)
        sheet = MonodomainSheet(cell, diff, dt=dt)

        def beat():
            for _ in range(round(stim_dur / dt)):
                cell.stimulate(left, 20.0)
                sheet.step()
            for _ in range(round((s1_bcl - stim_dur) / dt)):
                sheet.step()

        for _ in range(n_s1):
            beat()
        # wait DI then deliver S2, recording probe traces
        for _ in range(round(di / dt)):
            sheet.step()

        ts, v1s, v2s = [], [], []
        t = 0.0
        for _ in range(round(stim_dur / dt)):
            cell.stimulate(left, 20.0)
            sheet.step()
            t += dt
            ts.append(t); v1s.append(float(cell.V[0, p1])); v2s.append(float(cell.V[0, p2]))
        for _ in range(round(120.0 / dt)):
            sheet.step()
            t += dt
            ts.append(t); v1s.append(float(cell.V[0, p1])); v2s.append(float(cell.V[0, p2]))

        t1 = _upstroke_time(ts, v1s)
        t2 = _upstroke_time(ts, v2s)
        cv = dist_cm / (t2 - t1) * 1000.0  # cm per s
        cvs.append(cv)
    return dis, cvs
