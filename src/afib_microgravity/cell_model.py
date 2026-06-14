"""Common interface so AtrialSheet can drive either the AP or CRN kinetics.

A CellModel owns its own state arrays (the membrane variable plus any gating /
concentration variables). reaction_step advances the *local* kinetics one dt and
returns the membrane time-derivative contribution (dV/dt from reaction only); it
does NOT mutate the membrane variable V (the sheet owns the V update so it can add
diffusion). The sheet writes the diffused membrane field back via set_V.
"""
from __future__ import annotations

from abc import ABC, abstractmethod

import numpy as np

from .model import APParams  # existing dataclass


class CellModel(ABC):
    @property
    @abstractmethod
    def V(self) -> np.ndarray: ...

    @abstractmethod
    def set_V(self, V: np.ndarray) -> None: ...

    @abstractmethod
    def reaction_step(self, dt: float) -> np.ndarray:
        """Advance local kinetics by dt; return dV/dt_reaction (same shape). Must NOT mutate V."""

    @abstractmethod
    def stimulate(self, mask: np.ndarray, value: float) -> None: ...


class APCell(CellModel):
    """Aliev-Panfilov kinetics behind the CellModel interface (dimensionless)."""

    def __init__(self, shape=(200, 200), params: APParams | None = None):
        self.p = params or APParams()
        self.u = np.zeros(shape, dtype=float)
        self.v = np.zeros(shape, dtype=float)

    @property
    def V(self) -> np.ndarray:
        return self.u

    def set_V(self, V: np.ndarray) -> None:
        self.u = V

    def reaction_step(self, dt: float) -> np.ndarray:
        u, v, p = self.u, self.v, self.p
        i_ion = p.k * u * (u - p.a) * (u - 1.0) + u * v
        eps = p.eps0 + p.mu1 * v / (u + p.mu2)
        dv = eps * (-v - p.k * u * (u - p.a - 1.0))
        self.v = v + dt * dv
        return -i_ion  # reaction contribution to du/dt (diffusion added by sheet)

    def stimulate(self, mask: np.ndarray, value: float = 1.0) -> None:
        self.u[mask] = value
