"""In-silico model of atrial fibrillation under spaceflight (microgravity) remodelling.

A compact reaction-diffusion testbed for the hypothesis that microgravity-induced
atrial remodelling lowers the threshold for, and sustains, fibrillatory re-entry.
"""

from .model import APParams, AtrialSheet
from .metrics import (
    count_phase_singularities,
    dominant_frequency,
    phase_field,
    planar_conduction_velocity,
)
from .protocols import (
    induce_spiral,
    planar_s1,
    cross_field_s2,
    seed_broken_wavefront,
)
from .remodeling import make_condition, fibrosis_field

__version__ = "0.1.0"

__all__ = [
    "APParams",
    "AtrialSheet",
    "count_phase_singularities",
    "dominant_frequency",
    "phase_field",
    "planar_conduction_velocity",
    "induce_spiral",
    "planar_s1",
    "cross_field_s2",
    "seed_broken_wavefront",
    "make_condition",
    "fibrosis_field",
]
