"""Microgravity-induced atrial remodelling, mapped onto model parameters.

This module is the scientific heart of the project: it translates what spaceflight
physiology does to the atrium into concrete changes in the Aliev-Panfilov substrate.
Each mapping is a *hypothesis*, deliberately simple and individually toggleable so the
contribution of each mechanism can be isolated.

Three documented spaceflight changes, and how we represent them:

1. **Atrial dilation.**  Headward fluid redistribution in microgravity raises atrial
   filling and stretches the chamber.  A larger tissue area, in wavelengths, leaves
   more room for multiple coexisting wavelets.  -> larger grid.

2. **Structural remodelling / fibrosis.**  Deconditioning and altered loading promote
   patchy interstitial fibrosis, i.e. regions of poor cell-to-cell coupling that
   fragment wavefronts.  -> spatially random patches of reduced ``D``.

3. **Electrical / autonomic remodelling.**  Microgravity shifts sympatho-vagal balance
   and shortens atrial refractoriness, steepening restitution and shortening the
   wavelength -- classic pro-fibrillatory changes.  -> slightly faster recovery
   (larger ``eps0``).

NB: these are hypothesis-generating caricatures, not validated quantitative claims.
See the README's "Scientific honesty" section.
"""

from __future__ import annotations

import numpy as np

from .model import APParams


def baseline_params() -> APParams:
    """Healthy, ground-level atrial kinetics."""
    return APParams()


def microgravity_params(autonomic_severity: float = 1.0) -> APParams:
    """Electrical remodelling: shortened refractoriness via faster recovery.

    ``autonomic_severity`` in [0, ~2] scales the effect; 0 recovers baseline.
    """
    p = APParams()
    p.eps0 = p.eps0 * (1.0 + 1.5 * autonomic_severity)
    return p


def fibrosis_field(shape, density: float = 0.0, d_healthy: float = 1.0,
                   d_scar: float = 0.05, patch: int = 6, seed: int = 0):
    """Build a diffusion field with random low-coupling fibrotic patches.

    Parameters
    ----------
    density : float
        Fraction of patch sites converted to scar (0 = homogeneous tissue).
    patch : int
        Side length of each square fibrotic patch (cells).
    """
    ny, nx = shape
    field = np.full(shape, d_healthy, dtype=float)
    if density <= 0:
        return field
    rng = np.random.default_rng(seed)
    n_patches = int(density * (ny * nx) / (patch * patch))
    for _ in range(n_patches):
        i = rng.integers(0, max(1, ny - patch))
        j = rng.integers(0, max(1, nx - patch))
        field[i:i + patch, j:j + patch] = d_scar
    return field


def make_condition(name: str, base_shape=(220, 220)):
    """Return ``(shape, params, D_field)`` for a named condition.

    Conditions
    ----------
    ``"ground"``      : baseline atrium, homogeneous, nominal size.
    ``"microgravity"``: dilated (larger), fibrotic, and electrically remodelled.
    """
    if name == "ground":
        shape = base_shape
        return shape, baseline_params(), fibrosis_field(shape, density=0.0)
    if name == "microgravity":
        # dilation: ~30% larger linear dimension
        shape = (int(base_shape[0] * 1.3), int(base_shape[1] * 1.3))
        params = microgravity_params(autonomic_severity=1.0)
        d_field = fibrosis_field(shape, density=0.55, seed=7)
        return shape, params, d_field
    raise ValueError(f"unknown condition {name!r}")


# ---------------------------------------------------------------------------
# CRN-based remodelling functions
# ---------------------------------------------------------------------------

import dataclasses

from .crn import CRNCell, CRNParams
from .diffusion import AnisotropicDiffusion
from .fibrosis import correlated_fibrosis


def baseline_crn_params() -> CRNParams:
    """Default CRN 1998 parameters representing healthy ground-level atrium."""
    return CRNParams()


def microgravity_crn_params(severity: float = 1.0) -> CRNParams:
    """AF-type electrical remodelling as observed in chronic atrial fibrillation
    and hypothesised post-spaceflight.

    Maps the substrate onto CRN ionic conductances following the canonical
    AF electrical-remodelling pattern (Courtemanche et al. 1999 "AF" variant;
    Bosch 1999 / Workman 2001 human-atrial AF voltage-clamp data), which is
    the change set needed to collapse the action-potential duration enough for
    re-entry to fit a feasible sheet:

    * **ICaL** (g_CaL) strongly reduced → loss of the plateau, the dominant
      APD-shortening lever (the canonical AF reduction is ~65-70%).
    * **Ito**  (g_to)  reduced → less early-repolarisation reserve.
    * **IKur** (g_Kur_scale) reduced → less sustained outward K current.
    * **IK1**  (g_K1)  increased → faster terminal repolarisation and a more
      hyperpolarised, excitable resting state (a hallmark of AF remodelling).
    * **INa**  (g_Na)  unchanged → CV preserved, so the wavelength shortens
      almost entirely through APD.

    Together these collapse atrial APD90 from ~300 ms (baseline) to ~120-160 ms,
    shortening the wavelength to a few centimetres so a rotor can be hosted on a
    laptop-feasible sheet — the electrophysiological signature of an
    AF-vulnerable substrate.

    Parameters
    ----------
    severity:
        0 → recovers baseline exactly; 1 → full canonical AF remodelling
        (ICaL -70%, Ito -50%, IKur -50%, IK1 +100%).
    """
    base = CRNParams()
    return dataclasses.replace(
        base,
        g_CaL=base.g_CaL * (1.0 - 0.70 * severity),
        g_to=base.g_to * (1.0 - 0.50 * severity),
        g_Kur_scale=base.g_Kur_scale * (1.0 - 0.50 * severity),
        g_K1=base.g_K1 * (1.0 + 1.0 * severity),
    )


def make_condition_crn(name: str, base_shape=(200, 200), seed: int = 0):
    """Return ``(cell, diffusion)`` for a named CRN-based condition.

    Conditions
    ----------
    ``"ground"``       : baseline CRN kinetics, homogeneous diffusion.
    ``"microgravity"`` : dilated tissue, AF-remodelled kinetics, fibrotic coupling.
    """
    # d_long calibrated so the planar conduction velocity is ~58 cm/s
    # (physiological human atrium, target band 55-65); d_trans gives a ~3:1
    # anisotropy ratio. See experiments/calibrate_cv.py for the CV(d_long) sweep.
    D_LONG = 0.15
    D_TRANS = 0.05
    if name == "ground":
        shape = base_shape
        cell = CRNCell(shape=shape, params=baseline_crn_params())
        diff = AnisotropicDiffusion(
            shape=shape,
            d_long=D_LONG,
            d_trans=D_TRANS,
            theta=np.zeros(shape),
            dx=0.25,
            coupling=None,
        )
        return cell, diff

    if name == "microgravity":
        shape = (int(base_shape[0] * 1.3), int(base_shape[1] * 1.3))
        cell = CRNCell(shape=shape, params=microgravity_crn_params(severity=1.0))
        coupling = correlated_fibrosis(shape, density=0.3, corr_len=4.0, seed=seed)
        diff = AnisotropicDiffusion(
            shape=shape,
            d_long=D_LONG,
            d_trans=D_TRANS,
            theta=np.zeros(shape),
            dx=0.25,
            coupling=coupling,
        )
        return cell, diff

    raise ValueError(f"unknown condition {name!r}")
