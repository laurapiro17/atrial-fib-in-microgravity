"""Anisotropic monodomain diffusion: div(D grad V) with a per-cell conductivity
tensor built from a fibre-angle field. No-flux (Neumann) boundaries.

D = R(theta) diag(d_long, d_trans) R(theta)^T, so
  Dxx = d_long cos^2 + d_trans sin^2
  Dyy = d_long sin^2 + d_trans cos^2
  Dxy = (d_long - d_trans) sin cos
The cross (Dxy) term uses centred differences; axial terms use a harmonic-mean
finite-volume flux (robust across sharp conductivity contrasts, e.g. fibrosis).
"""
from __future__ import annotations

import numpy as np


class AnisotropicDiffusion:
    def __init__(self, shape, d_long, d_trans, theta, dx=0.25, coupling=None):
        self.ny, self.nx = shape
        self.dx = float(dx)
        c, s = np.cos(theta), np.sin(theta)
        scale = np.ones(shape) if coupling is None else np.asarray(coupling, float)
        dl = d_long * scale
        dt_ = d_trans * scale
        self.Dxx = dl * c * c + dt_ * s * s
        self.Dyy = dl * s * s + dt_ * c * c
        self.Dxy = (dl - dt_) * s * c

    @staticmethod
    def _hmean(a, b):
        s = a + b
        out = np.zeros_like(a)
        nz = s > 0
        out[nz] = 2.0 * a[nz] * b[nz] / s[nz]
        return out

    def apply(self, V):
        inv = 1.0 / (self.dx * self.dx)
        Dr = self._hmean(self.Dxx, np.roll(self.Dxx, -1, 1))
        fx_r = Dr * (np.roll(V, -1, 1) - V); fx_r[:, -1] = 0.0
        Dl = self._hmean(self.Dxx, np.roll(self.Dxx, 1, 1))
        fx_l = Dl * (V - np.roll(V, 1, 1)); fx_l[:, 0] = 0.0
        Dd = self._hmean(self.Dyy, np.roll(self.Dyy, -1, 0))
        fy_d = Dd * (np.roll(V, -1, 0) - V); fy_d[-1, :] = 0.0
        Du = self._hmean(self.Dyy, np.roll(self.Dyy, 1, 0))
        fy_u = Du * (V - np.roll(V, 1, 0)); fy_u[0, :] = 0.0
        axial = ((fx_r - fx_l) + (fy_d - fy_u)) * inv
        dVdy = (np.roll(V, -1, 0) - np.roll(V, 1, 0)) * 0.5
        cross = (np.roll(self.Dxy * dVdy, -1, 1) -
                 np.roll(self.Dxy * dVdy, 1, 1)) * 0.5 * inv
        cross[0, :] = cross[-1, :] = cross[:, 0] = cross[:, -1] = 0.0
        return axial + cross

    @property
    def Dmax(self):
        return float(max(self.Dxx.max(), self.Dyy.max()))
