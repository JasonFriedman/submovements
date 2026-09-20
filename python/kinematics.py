"""Shared kinematic preparation for plotting decomposed submovements."""

from __future__ import annotations

import numpy as np

from min_jerk import (
    minimum_jerk_position_1d,
    minimum_jerk_position_2d,
    minimum_jerk_position_3d,
    minimum_jerk_velocity_1d,
    minimum_jerk_velocity_2d,
    minimum_jerk_velocity_3d,
)


def prepare_submovement_kinematics(parameters, dimensions: int, t=None, initial_positions=None):
    if dimensions not in (1, 2, 3):
        raise ValueError('dimensions must be 1, 2, or 3')

    parameters = np.asarray(parameters, dtype=float).ravel()
    pps = dimensions + 2
    if parameters.size % pps != 0:
        raise ValueError(f'The parameters vector must have a length that is a multiple of {pps}')

    num = parameters.size // pps
    matrix = parameters.reshape(num, pps)
    matrix = matrix[np.argsort(matrix[:, 0]), :]

    t0 = matrix[:, 0]
    duration = matrix[:, 1]
    amp = matrix[:, 2:]

    if initial_positions is None:
        initial_positions = np.zeros(dimensions)
    initial_positions = np.asarray(initial_positions, dtype=float).reshape(dimensions)

    starts = np.zeros((num, dimensions), dtype=float)
    starts[0, :] = initial_positions
    if num > 1:
        for d in range(dimensions):
            starts[1:, d] = starts[0, d] + np.cumsum(amp[:-1, d])

    if t is None:
        t = np.linspace(np.min(t0), np.max(t0 + duration), 100)
    t = np.asarray(t, dtype=float)

    vel = np.zeros((num, t.size, dimensions), dtype=float)
    pos = np.zeros((num, t.size, dimensions), dtype=float)

    for s in range(num):
        if dimensions == 1:
            vel[s, :, 0] = minimum_jerk_velocity_1d(t0[s], duration[s], amp[s, 0], t)
            pos[s, :, 0] = minimum_jerk_position_1d(t0[s], duration[s], amp[s, 0], starts[s, 0], t)
        elif dimensions == 2:
            vx, vy, _ = minimum_jerk_velocity_2d(t0[s], duration[s], amp[s, 0], amp[s, 1], t)
            x, y = minimum_jerk_position_2d(t0[s], duration[s], amp[s, 0], amp[s, 1], starts[s, 0], starts[s, 1], t)
            vel[s, :, 0], vel[s, :, 1] = vx, vy
            pos[s, :, 0], pos[s, :, 1] = x, y
        else:
            vx, vy, vz, _ = minimum_jerk_velocity_3d(t0[s], duration[s], amp[s, 0], amp[s, 1], amp[s, 2], t)
            x, y, z = minimum_jerk_position_3d(
                t0[s], duration[s], amp[s, 0], amp[s, 1], amp[s, 2], starts[s, 0], starts[s, 1], starts[s, 2], t
            )
            vel[s, :, 0], vel[s, :, 1], vel[s, :, 2] = vx, vy, vz
            pos[s, :, 0], pos[s, :, 1], pos[s, :, 2] = x, y, z

    pos_relative = pos.copy()
    for s in range(1, num):
        pos_relative[s, :, :] -= starts[s, :]

    sum_pos = np.sum(pos_relative, axis=0)

    return {
        'numSubmovements': num,
        't0': t0,
        'D': duration,
        'A': amp,
        'starts': starts,
        't': t,
        'vel': vel,
        'pos': pos,
        'sumPos': sum_pos,
    }
