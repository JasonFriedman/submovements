"""Minimum-jerk basis and kinematic utilities (1D/2D/3D)."""

from __future__ import annotations

import numpy as np


def minimum_jerk_basis(t0: float, duration: float, t: np.ndarray):
    t = np.asarray(t, dtype=float)
    nt = (t - t0) / duration
    active = (nt >= 0) & (nt <= 1)
    before = nt < 0
    after = nt > 1

    position_basis = -15 * nt**4 + 6 * nt**5 + 10 * nt**3
    velocity_basis = 30 * nt**2 - 60 * nt**3 + 30 * nt**4

    return nt, active, before, after, position_basis, velocity_basis


def minimum_jerk_position_1d(t0: float, duration: float, amplitude: float, x0: float, t: np.ndarray):
    nt, active, before, after, position_basis, _ = minimum_jerk_basis(t0, duration, t)
    x = np.zeros_like(nt)
    x[active] = amplitude * position_basis[active] + x0
    x[before] = x0
    x[after] = x0 + amplitude
    return x


def minimum_jerk_velocity_1d(t0: float, duration: float, amplitude: float, t: np.ndarray):
    nt, active, _, _, _, velocity_basis = minimum_jerk_basis(t0, duration, t)
    velocity = np.zeros_like(nt)
    velocity[active] = amplitude / duration * velocity_basis[active]
    return velocity


def minimum_jerk_acceleration_1d(t0: float, duration: float, amplitude: float, t: np.ndarray):
    nt, active, _, _, _, _ = minimum_jerk_basis(t0, duration, t)
    acc = np.zeros_like(nt)
    acc[active] = amplitude / (duration**2) * (60 * nt[active] - 180 * nt[active] ** 2 + 120 * nt[active] ** 3)
    return acc


def minimum_jerk_jerk_1d(t0: float, duration: float, amplitude: float, t: np.ndarray):
    nt, active, _, _, _, _ = minimum_jerk_basis(t0, duration, t)
    jerk = np.zeros_like(nt)
    jerk[active] = amplitude / (duration**3) * (60 - 360 * nt[active] + 360 * nt[active] ** 2)
    return jerk


def minimum_jerk_position_2d(t0: float, duration: float, ax: float, ay: float, x0: float, y0: float, t: np.ndarray):
    x = minimum_jerk_position_1d(t0, duration, ax, x0, t)
    y = minimum_jerk_position_1d(t0, duration, ay, y0, t)
    return x, y


def minimum_jerk_velocity_2d(t0: float, duration: float, ax: float, ay: float, t: np.ndarray):
    vx = minimum_jerk_velocity_1d(t0, duration, ax, t)
    vy = minimum_jerk_velocity_1d(t0, duration, ay, t)
    speed = minimum_jerk_velocity_1d(t0, duration, np.sqrt(ax**2 + ay**2), t)
    return vx, vy, speed


def minimum_jerk_position_3d(
    t0: float,
    duration: float,
    ax: float,
    ay: float,
    az: float,
    x0: float,
    y0: float,
    z0: float,
    t: np.ndarray,
):
    x = minimum_jerk_position_1d(t0, duration, ax, x0, t)
    y = minimum_jerk_position_1d(t0, duration, ay, y0, t)
    z = minimum_jerk_position_1d(t0, duration, az, z0, t)
    return x, y, z


def minimum_jerk_velocity_3d(t0: float, duration: float, ax: float, ay: float, az: float, t: np.ndarray):
    vx = minimum_jerk_velocity_1d(t0, duration, ax, t)
    vy = minimum_jerk_velocity_1d(t0, duration, ay, t)
    vz = minimum_jerk_velocity_1d(t0, duration, az, t)
    speed = minimum_jerk_velocity_1d(t0, duration, np.sqrt(ax**2 + ay**2 + az**2), t)
    return vx, vy, vz, speed
