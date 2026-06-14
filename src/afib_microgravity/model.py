"""Aliev-Panfilov monodomain model of a 2D atrial tissue sheet.

The Aliev-Panfilov (1996) model is a two-variable, dimensionless reaction-diffusion
system that reproduces the essential excitable dynamics of cardiac tissue --
including the spiral-wave re-entry that underlies atrial fibrillation -- at a tiny
fraction of the cost of an ionic model.  We use it here as the substrate on which to
ask: *how does the structural and electrical remodelling produced by microgravity
change re-entry dynamics?*

Governing equations (per unit dimensionless area):

    du/dt = div(D grad u) - k u (u - a)(u - 1) - u v
    dv/dt = (eps0 + mu1 v / (u + mu2)) (-v - k u (u - a - 1))

where ``u`` is the (fast) transmembrane potential and ``v`` a (slow) recovery
variable.  ``D`` is allowed to vary in space so that fibrosis can be represented as
patches of reduced cell-to-cell coupling.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class APParams:
    """Aliev-Panfilov kinetic parameters.

    The defaults are the canonical values from Aliev & Panfilov (1996) and reliably
    produce sustained spiral waves on a sheet a few wavelengths across.
    """

    a: float = 0.05      # excitation threshold
    k: float = 8.0       # excitation strength
    eps0: float = 0.002  # base recovery rate
    mu1: float = 0.2     # recovery coupling
    mu2: float = 0.3     # recovery coupling
    D: float = 1.0       # default (uniform) diffusion coefficient


class AtrialSheet:
    """A 2D excitable sheet integrated with explicit Euler in time.

    Parameters
    ----------
    shape : (ny, nx)
        Grid size in cells.
    params : APParams
        Kinetic parameters.
    D_field : ndarray, optional
        Per-cell diffusion coefficient.  Defaults to a uniform ``params.D``.
        Spatial heterogeneity here is how we model fibrosis.
    dx, dt : float
        Space and time steps (dimensionless).  ``dt`` must satisfy the explicit
        diffusion-stability bound ``dt <= dx**2 / (4 * Dmax)``.
    """

    def __init__(self, shape=(200, 200), params: APParams | None = None,
                 D_field=None, dx: float = 1.0, dt: float = 0.02, seed: int = 0):
        self.ny, self.nx = shape
        self.p = params or APParams()
        self.dx = float(dx)
        self.dt = float(dt)
        self.u = np.zeros(shape, dtype=float)
        self.v = np.zeros(shape, dtype=float)
        if D_field is None:
            self.D_field = np.full(shape, self.p.D, dtype=float)
        else:
            self.D_field = np.asarray(D_field, dtype=float)
            if self.D_field.shape != shape:
                raise ValueError("D_field shape must match grid shape")
        self.t = 0.0
        self._check_stability()

    def _check_stability(self):
        bound = self.dx * self.dx / (4.0 * max(self.D_field.max(), 1e-12))
        if self.dt > bound:
            raise ValueError(
                f"dt={self.dt} violates diffusion stability bound {bound:.4g}; "
                "reduce dt or D."
            )

    @staticmethod
    def _hmean(a, b):
        """Harmonic mean of neighbouring diffusivities (good across sharp contrasts)."""
        s = a + b
        out = np.zeros_like(a)
        nz = s > 0
        out[nz] = 2.0 * a[nz] * b[nz] / s[nz]
        return out

    def _div_D_grad(self, u):
        """Finite-volume div(D grad u) with no-flux (Neumann) boundaries.

        Face diffusivities use the harmonic mean of the two adjacent cells, and the
        flux across every domain boundary face is set to zero.  On a uniform field
        this reduces to the standard 5-point Laplacian.
        """
        D = self.D_field
        inv_dx2 = 1.0 / (self.dx * self.dx)

        # +x faces (between column i and i+1)
        Dr = self._hmean(D, np.roll(D, -1, axis=1))
        fx_r = Dr * (np.roll(u, -1, axis=1) - u)
        fx_r[:, -1] = 0.0
        # -x faces
        Dl = self._hmean(D, np.roll(D, 1, axis=1))
        fx_l = Dl * (u - np.roll(u, 1, axis=1))
        fx_l[:, 0] = 0.0

        # +y faces
        Dd = self._hmean(D, np.roll(D, -1, axis=0))
        fy_d = Dd * (np.roll(u, -1, axis=0) - u)
        fy_d[-1, :] = 0.0
        # -y faces
        Du = self._hmean(D, np.roll(D, 1, axis=0))
        fy_u = Du * (u - np.roll(u, 1, axis=0))
        fy_u[0, :] = 0.0

        return ((fx_r - fx_l) + (fy_d - fy_u)) * inv_dx2

    def step(self):
        """Advance one explicit-Euler time step."""
        u, v, p = self.u, self.v, self.p
        diffusion = self._div_D_grad(u)
        i_ion = p.k * u * (u - p.a) * (u - 1.0) + u * v
        eps = p.eps0 + p.mu1 * v / (u + p.mu2)
        du = diffusion - i_ion
        dv = eps * (-v - p.k * u * (u - p.a - 1.0))
        self.u = u + self.dt * du
        self.v = v + self.dt * dv
        self.t += self.dt

    def run(self, n_steps: int):
        """Advance ``n_steps`` steps in place."""
        for _ in range(n_steps):
            self.step()

    def stimulate(self, mask, value: float = 1.0):
        """Raise ``u`` to ``value`` wherever ``mask`` is True (a stimulus)."""
        self.u[mask] = value


class MonodomainSheet:
    """Operator-split monodomain sheet: explicit reaction (CellModel) then explicit
    anisotropic diffusion. The cell model returns dV/dt from reaction and does not
    mutate V; this sheet owns the membrane update V += dt*(reaction + diffusion).
    Stimulate via cell.stimulate(mask, value) before/between steps."""

    def __init__(self, cell, diffusion, dt=0.01):
        self.cell = cell
        self.diff = diffusion
        self.dt = float(dt)
        self.t = 0.0
        self._check_stability()

    def _check_stability(self):
        bound = self.diff.dx ** 2 / (4.0 * max(self.diff.Dmax, 1e-12))
        if self.dt > bound:
            raise ValueError(
                f"dt={self.dt} exceeds explicit diffusion stability bound {bound:.4g}"
            )

    def step(self):
        dV_react = self.cell.reaction_step(self.dt)
        V = self.cell.V + self.dt * dV_react
        V = V + self.dt * self.diff.apply(V)
        self.cell.set_V(V)
        self.t += self.dt

    def run(self, n_steps):
        for _ in range(int(n_steps)):
            self.step()
