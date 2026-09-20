"""Compatibility API for movement decomposition and plotting.

This module keeps the historic function names used by notebooks while delegating
implementation to the refactored modular Python files.
"""

from __future__ import annotations

import numpy as np

from decompose import decompose_1d, decompose_2d, decompose_3d
from plotting import (
    plot_position,
    plot_velocity,
    plot_velocity_and_reconstruction_2d,
    plot_submovements_1d,
    plot_submovements_2d,
    plot_submovements_3d,
)
from windows import decompose_1d_windows, decompose_2d_windows, decompose_3d_windows


def decompose_1D(time: np.ndarray, vel: np.ndarray, n_sub_movement: int = 4, a_rng: tuple[float, float] = (-5.0, 5.0)):
    """Legacy wrapper for 1D decomposition."""
    best_error, best_parameters, _ = decompose_1d(time, vel, num_submovements=n_sub_movement, a_rng=a_rng)
    return best_error, np.asarray(best_parameters).reshape(-1, 3)


def decompose_2D(
    time: np.ndarray,
    vel: np.ndarray,
    n_sub_movement: int = 4,
    x_rng: tuple[float, float] = (-5.0, 5.0),
    y_rng: tuple[float, float] = (0.1, 5.0),
):
    """Legacy wrapper for 2D decomposition."""
    best_error, best_parameters, _ = decompose_2d(
        time,
        vel,
        num_submovements=n_sub_movement,
        x_rng=x_rng,
        y_rng=y_rng,
    )
    return best_error, np.asarray(best_parameters).reshape(-1, 4)


def decompose_3D(
    time: np.ndarray,
    vel: np.ndarray,
    n_sub_movement: int = 4,
    x_rng: tuple[float, float] = (-5.0, 5.0),
    y_rng: tuple[float, float] = (0.1, 5.0),
    z_rng: tuple[float, float] = (-5.0, 5.0),
):
    """Legacy wrapper for 3D decomposition."""
    best_error, best_parameters, _ = decompose_3d(
        time,
        vel,
        num_submovements=n_sub_movement,
        x_rng=x_rng,
        y_rng=y_rng,
        z_rng=z_rng,
    )
    return best_error, np.asarray(best_parameters).reshape(-1, 5)


def plot_submovements_1D(
    parameters, 
    t: np.ndarray = None, 
    plot_type: int = 1, 
    x0: float = 0.0
):
    return plot_submovements_1d(parameters, t=t, plot_type=plot_type, x0=x0)


def plot_submovements_2D(
    parameters, 
    t: np.ndarray = None, 
    plot_type: int = 1, 
    x0: float = 0.0, 
    y0: float = 0.0
):
    return plot_submovements_2d(parameters, t=t, plot_type=plot_type, x0=x0, y0=y0)


def plot_velocity_and_reconstruction_2D(
    velocity: np.ndarray,
    parameters,
    t: np.ndarray = None,
    x0: float = 0.0,
    y0: float = 0.0,
):
    return plot_velocity_and_reconstruction_2d(velocity, parameters, t=t, x0=x0, y0=y0)


def plot_submovements_3D(
    parameters,
    t: np.ndarray = None,
    plot_type: int = 1,
    x0: float = 0.0,
    y0: float = 0.0,
    z0: float = 0.0,
):
    return plot_submovements_3d(parameters, t=t, plot_type=plot_type, x0=x0, y0=y0, z0=z0)


__all__ = [
    'plot_position',
    'plot_velocity',
    'decompose_1D',
    'decompose_2D',
    'decompose_3D',
    'plot_submovements_1D',
    'plot_submovements_2D',
    'plot_submovements_3D',
    'plot_velocity_and_reconstruction_2D',
    'decompose_1d',
    'decompose_2d',
    'decompose_3d',
    'decompose_1d_windows',
    'decompose_2d_windows',
    'decompose_3d_windows',
]
