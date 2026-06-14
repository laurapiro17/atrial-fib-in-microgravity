"""Gravitational scaling law of atrial arrhythmogenesis.

A continuous map from gravitational level ``g`` (in Earth-g units: g=1.0 is
Earth, g=0 is free fall) to the microgravity remodelling knobs, plus the
dimensionless Cardiac Gravitational Number

    N_g = L(g) / WL(g) = L0 * dilation(g) / (CV * APD90(g) / 1000)

Re-entry can be hosted when N_g >= 1. The critical gravity g* solves
N_g(g*) = 1.

HONESTY: the g -> remodelling map is a deliberately simple *linear* hypothesis
in the cephalad fluid-shift drive (1 - g), calibrated so g=1 recovers the
validated ground baseline and g=0 recovers the project's existing microgravity
operating point (see remodeling.py). The linear form and the mechano-electric
sensitivities are assumptions swept in sensitivity analysis, not measured
facts. Evidence base: Khine 2018 (10.1161/CIRCEP.117.005959, atrial
enlargement + AF risk markers in spaceflight) and Ravelli 2003
(10.1016/s0079-6107(03)00011-7, stretch shortens refractoriness / slows CV).
"""
from __future__ import annotations

from dataclasses import dataclass

EARTH_G = 1.0
# Endpoints recovered at g=0, taken from remodeling.make_condition_crn("microgravity").
SEVERITY_MAX = 1.0
DILATION_MAX = 1.3
DENSITY_MAX = 0.3


@dataclass(frozen=True)
class GravityRemodeling:
    gravity: float
    severity: float
    dilation: float
    density: float


def gravity_to_remodeling(g: float) -> GravityRemodeling:
    """Map gravitational level ``g`` to remodelling knobs.

    Linear in the fluid-shift drive ``(1 - g)``; ``g`` clamped to [0, 1].
    """
    g = min(max(g, 0.0), 1.0)
    drive = 1.0 - g
    return GravityRemodeling(
        gravity=g,
        severity=SEVERITY_MAX * drive,
        dilation=1.0 + (DILATION_MAX - 1.0) * drive,
        density=DENSITY_MAX * drive,
    )


def wavelength_cm(apd90_ms: float, cv_cm_s: float = 58.0) -> float:
    """Electrical wavelength WL = CV * APD90, in cm."""
    return cv_cm_s * (apd90_ms / 1000.0)


def cardiac_gravitational_number(
    apd90_ms: float,
    dilation: float,
    l0_cm: float = 8.0,
    cv_cm_s: float = 58.0,
) -> float:
    """Dimensionless N_g = L0*dilation / WL.

    l0_cm: characteristic atrial path length at ground (human left-atrial
    circumference ~8-12 cm; default 8). It dilates with the fluid shift.
    """
    wl = wavelength_cm(apd90_ms, cv_cm_s)
    return (l0_cm * dilation) / wl


def interpolate_crossing(xs, ys, target: float = 1.0):
    """First x (linear-interpolated) where y crosses ``target``.

    ``xs`` ascending. Returns None if no sign change of (y - target).
    """
    for i in range(len(xs) - 1):
        d0 = ys[i] - target
        d1 = ys[i + 1] - target
        if d0 == 0.0:
            return xs[i]
        if d0 * d1 < 0.0:
            t = (target - ys[i]) / (ys[i + 1] - ys[i])
            return xs[i] + t * (xs[i + 1] - xs[i])
    return None
