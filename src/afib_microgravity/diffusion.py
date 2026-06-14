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

from .crn_numba import HAVE_NUMBA

if HAVE_NUMBA:
    from numba import njit, prange

    @njit(cache=True, parallel=True, fastmath=False)
    def _apply_kernel(V, Dxx, Dyy, Dxy, inv):
        ny, nx = V.shape
        out = np.zeros((ny, nx), dtype=V.dtype)
        for i in prange(ny):
            for jx in range(nx):
                Vc = V[i, jx]
                Dc_x = Dxx[i, jx]
                Dc_y = Dyy[i, jx]
                acc = 0.0
                # +x face (between jx and jx+1)
                if jx + 1 < nx:
                    Dr = Dxx[i, jx + 1]
                    s = Dc_x + Dr
                    if s > 0.0:
                        acc += (2.0 * Dc_x * Dr / s) * (V[i, jx + 1] - Vc)
                # -x face
                if jx - 1 >= 0:
                    Dl = Dxx[i, jx - 1]
                    s = Dc_x + Dl
                    if s > 0.0:
                        acc -= (2.0 * Dc_x * Dl / s) * (Vc - V[i, jx - 1])
                # +y face
                if i + 1 < ny:
                    Dd = Dyy[i + 1, jx]
                    s = Dc_y + Dd
                    if s > 0.0:
                        acc += (2.0 * Dc_y * Dd / s) * (V[i + 1, jx] - Vc)
                # -y face
                if i - 1 >= 0:
                    Du = Dyy[i - 1, jx]
                    s = Dc_y + Du
                    if s > 0.0:
                        acc -= (2.0 * Dc_y * Du / s) * (Vc - V[i - 1, jx])
                out[i, jx] = acc * inv
        # cross (Dxy) term: matches the NumPy centred-difference scheme, with the
        # one-cell-wide border zeroed exactly as in the reference implementation.
        for i in prange(1, ny - 1):
            for jx in range(1, nx - 1):
                dVdy_r = 0.5 * (V[i + 1, jx + 1] - V[i - 1, jx + 1])
                dVdy_l = 0.5 * (V[i + 1, jx - 1] - V[i - 1, jx - 1])
                out[i, jx] += (0.5 * inv
                               * (Dxy[i, jx + 1] * dVdy_r
                                  - Dxy[i, jx - 1] * dVdy_l))
        return out


class AnisotropicDiffusion:
    def __init__(self, shape, d_long, d_trans, theta, dx=0.25, coupling=None,
                 use_numba=None):
        self.ny, self.nx = shape
        self.dx = float(dx)
        self.use_numba = HAVE_NUMBA if use_numba is None else bool(use_numba)
        c, s = np.cos(theta), np.sin(theta)
        scale = np.ones(shape) if coupling is None else np.asarray(coupling, float)
        dl = d_long * scale
        dt_ = d_trans * scale
        self.Dxx = np.ascontiguousarray(dl * c * c + dt_ * s * s)
        self.Dyy = np.ascontiguousarray(dl * s * s + dt_ * c * c)
        self.Dxy = np.ascontiguousarray((dl - dt_) * s * c)

    @staticmethod
    def _hmean(a, b):
        s = a + b
        out = np.zeros_like(a)
        nz = s > 0
        out[nz] = 2.0 * a[nz] * b[nz] / s[nz]
        return out

    def apply(self, V):
        inv = 1.0 / (self.dx * self.dx)
        if self.use_numba:
            return _apply_kernel(np.ascontiguousarray(V),
                                 self.Dxx, self.Dyy, self.Dxy, inv)
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
