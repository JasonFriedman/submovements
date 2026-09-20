"""Decomposition APIs for 1D/2D/3D minimum-jerk models."""

from __future__ import annotations

import numpy as np

from constraints import remove_non_optimizer_constraints, resolve_fitting_constraints, resolve_pair_cancellation_penalty
from decomposition_engine import decompose_nd
from min_jerk import minimum_jerk_velocity_1d, minimum_jerk_velocity_2d, minimum_jerk_velocity_3d


def _pair_cancellation_penalty(parameters: np.ndarray, parameters_per_sub: int, amplitude_dimensions: int, settings: dict) -> float:
    weight = settings['weight']
    if weight <= 0:
        return 0.0

    matrix = parameters.reshape(-1, parameters_per_sub)
    t0 = matrix[:, 0]
    duration = matrix[:, 1]
    amplitude = matrix[:, 2 : 2 + amplitude_dimensions]

    penalty = 0.0
    for i in range(matrix.shape[0] - 1):
        for j in range(i + 1, matrix.shape[0]):
            dt = t0[i] - t0[j]
            dd = duration[i] - duration[j]
            da = amplitude[i] - amplitude[j]
            similarity = np.exp(-((dt / settings['t0Sigma']) ** 2 + (dd / settings['dSigma']) ** 2))
            penalty += weight * similarity * float(np.dot(da, da))

    return penalty


def _extend_if_necessary(parameters: np.ndarray, pps: int, time: np.ndarray, vel: np.ndarray, timedelta: float):
    matrix = parameters.reshape(-1, pps)
    last_time = float(np.max(matrix[:, 0] + matrix[:, 1]))
    last_time = np.round(last_time * (1 / timedelta)) / (1 / timedelta)
    if last_time <= time[-1]:
        return time, vel

    extension = np.arange(time[-1], last_time + timedelta, timedelta)
    extended_time = np.concatenate([time[:-1], extension])
    pad_length = extended_time.size - vel.shape[0]
    if vel.ndim == 1:
        vel = np.concatenate([vel, np.zeros(pad_length)])
    else:
        vel = np.vstack([vel, np.zeros((pad_length, vel.shape[1]))])
    return extended_time, vel


def calculate_error_mj_1d(parameters: np.ndarray, time: np.ndarray, vel: np.ndarray, timedelta: float = 0.005, penalty_settings: dict | None = None):
    settings = {'weight': 0.0, 't0Sigma': 0.05, 'dSigma': 0.05} if penalty_settings is None else penalty_settings

    time_eval, vel_eval = _extend_if_necessary(parameters, 3, time, vel, timedelta)
    trajectory = vel_eval[:, 0] if vel_eval.ndim == 2 else vel_eval

    matrix = parameters.reshape(-1, 3)
    predicted = np.zeros((matrix.shape[0], time_eval.size))
    for i, (t0, duration, amplitude) in enumerate(matrix):
        active = (time_eval >= t0) & (time_eval <= t0 + duration)
        predicted[i, active] = minimum_jerk_velocity_1d(t0, duration, amplitude, time_eval[active])

    sum_predicted = np.sum(predicted, axis=0)
    sum_traj_sq = np.sum(trajectory**2)
    if sum_traj_sq == 0:
        sum_traj_sq = 1.0

    epsilon = np.sum((sum_predicted - trajectory) ** 2) / sum_traj_sq
    epsilon += _pair_cancellation_penalty(parameters, 3, 1, settings)

    return epsilon, sum_predicted


def calculate_error_mj_2d(parameters: np.ndarray, time: np.ndarray, vel: np.ndarray, timedelta: float = 0.005, penalty_settings: dict | None = None):
    settings = {'weight': 0.0, 't0Sigma': 0.05, 'dSigma': 0.05} if penalty_settings is None else penalty_settings

    time_eval, vel_eval = _extend_if_necessary(parameters, 4, time, vel, timedelta)
    trajectory_x = vel_eval[:, 0]
    trajectory_y = vel_eval[:, 1]
    tangential = np.sqrt(trajectory_x**2 + trajectory_y**2)

    matrix = parameters.reshape(-1, 4)
    px = np.zeros((matrix.shape[0], time_eval.size))
    py = np.zeros((matrix.shape[0], time_eval.size))
    pt = np.zeros((matrix.shape[0], time_eval.size))

    for i, (t0, duration, ax, ay) in enumerate(matrix):
        active = (time_eval >= t0) & (time_eval <= t0 + duration)
        vx, vy, vt = minimum_jerk_velocity_2d(t0, duration, ax, ay, time_eval[active])
        px[i, active] = vx
        py[i, active] = vy
        pt[i, active] = vt

    sumx = np.sum(px, axis=0)
    sumy = np.sum(py, axis=0)
    sumt = np.sum(pt, axis=0)

    sum_traj_sq = np.sum(trajectory_x**2 + trajectory_y**2 + tangential**2)
    if sum_traj_sq == 0:
        sum_traj_sq = 1.0

    epsilon = np.sum((sumx - trajectory_x) ** 2 + (sumy - trajectory_y) ** 2 + (sumt - tangential) ** 2) / sum_traj_sq
    epsilon += _pair_cancellation_penalty(parameters, 4, 2, settings)

    return epsilon, np.column_stack([sumx, sumy])


def calculate_error_mj_3d(parameters: np.ndarray, time: np.ndarray, vel: np.ndarray, timedelta: float = 0.005, penalty_settings: dict | None = None):
    settings = {'weight': 0.0, 't0Sigma': 0.05, 'dSigma': 0.05} if penalty_settings is None else penalty_settings

    time_eval, vel_eval = _extend_if_necessary(parameters, 5, time, vel, timedelta)
    trajectory_x = vel_eval[:, 0]
    trajectory_y = vel_eval[:, 1]
    trajectory_z = vel_eval[:, 2]
    tangential = np.sqrt(trajectory_x**2 + trajectory_y**2 + trajectory_z**2)

    matrix = parameters.reshape(-1, 5)
    px = np.zeros((matrix.shape[0], time_eval.size))
    py = np.zeros((matrix.shape[0], time_eval.size))
    pz = np.zeros((matrix.shape[0], time_eval.size))
    pt = np.zeros((matrix.shape[0], time_eval.size))

    for i, (t0, duration, ax, ay, az) in enumerate(matrix):
        active = (time_eval >= t0) & (time_eval <= t0 + duration)
        vx, vy, vz, vt = minimum_jerk_velocity_3d(t0, duration, ax, ay, az, time_eval[active])
        px[i, active] = vx
        py[i, active] = vy
        pz[i, active] = vz
        pt[i, active] = vt

    sumx = np.sum(px, axis=0)
    sumy = np.sum(py, axis=0)
    sumz = np.sum(pz, axis=0)
    sumt = np.sum(pt, axis=0)

    sum_traj_sq = np.sum(trajectory_x**2 + trajectory_y**2 + trajectory_z**2 + tangential**2)
    if sum_traj_sq == 0:
        sum_traj_sq = 1.0

    epsilon = np.sum(
        (sumx - trajectory_x) ** 2
        + (sumy - trajectory_y) ** 2
        + (sumz - trajectory_z) ** 2
        + (sumt - tangential) ** 2
    ) / sum_traj_sq
    epsilon += _pair_cancellation_penalty(parameters, 5, 3, settings)

    return epsilon, np.column_stack([sumx, sumy, sumz])


def _decompose_generic(time, vel, num_submovements, lb0, ub0, pps, error_fn, fit_fn, fitting_constraints):
    if np.asarray(time).ndim != 1:
        raise ValueError('time must be 1D')
    if len(time) != len(vel):
        raise ValueError('time and vel must have equal length')

    if num_submovements is None:
        num_submovements = [1, 2, 3, 4]

    if isinstance(num_submovements, (list, tuple, np.ndarray)):
        if len(num_submovements) == 1:
            num_submovements = int(num_submovements[0])
        else:
            errors, parameters, fits = [], [], []
            for count in num_submovements:
                e, p, f = _decompose_generic(time, vel, int(count), lb0, ub0, pps, error_fn, fit_fn, fitting_constraints)
                errors.append(e)
                parameters.append(p)
                fits.append(f)
            return np.array(errors), parameters, fits

    return decompose_nd(
        np.asarray(time),
        int(num_submovements),
        np.asarray(lb0, dtype=float),
        np.asarray(ub0, dtype=float),
        pps,
        error_fn,
        fit_fn,
        fitting_constraints,
    )


def decompose_1d(time, vel, num_submovements=None, a_rng=(-5.0, 5.0), criteria=-np.inf, fitting_constraints=None):
    fitting_constraints = {} if fitting_constraints is None else dict(fitting_constraints)
    optimizer_constraints = remove_non_optimizer_constraints(fitting_constraints)
    constraints = resolve_fitting_constraints(optimizer_constraints)
    penalty = resolve_pair_cancellation_penalty(fitting_constraints)

    time = np.asarray(time, dtype=float)
    vel = np.asarray(vel, dtype=float)
    if vel.ndim == 2 and vel.shape[1] == 1:
        vel = vel[:, 0]

    lb0 = [0, constraints['minDuration'], a_rng[0]]
    ub0 = [max(time[-1] - constraints['minOnsetSpacing'], constraints['minUpperBoundTime']), constraints['maxDuration'], a_rng[1]]
    dt = float(time[1] - time[0]) if time.size > 1 else 0.005

    def err_fn(parameters):
        return calculate_error_mj_1d(parameters, time, vel, dt, penalty)[0]

    def fit_fn(parameters):
        return calculate_error_mj_1d(parameters, time, vel, dt, penalty)[1]

    return _decompose_generic(time, vel, num_submovements, lb0, ub0, 3, err_fn, fit_fn, optimizer_constraints)


def decompose_2d(time, vel, num_submovements=None, x_rng=(-5.0, 5.0), y_rng=(0.1, 5.0), criteria=-np.inf, fitting_constraints=None):
    fitting_constraints = {} if fitting_constraints is None else dict(fitting_constraints)
    optimizer_constraints = remove_non_optimizer_constraints(fitting_constraints)
    constraints = resolve_fitting_constraints(optimizer_constraints)
    penalty = resolve_pair_cancellation_penalty(fitting_constraints)

    time = np.asarray(time, dtype=float)
    vel = np.asarray(vel, dtype=float)

    lb0 = [0, constraints['minDuration'], x_rng[0], y_rng[0]]
    ub0 = [max(time[-1] - constraints['minOnsetSpacing'], constraints['minUpperBoundTime']), constraints['maxDuration'], x_rng[1], y_rng[1]]
    dt = float(time[1] - time[0]) if time.size > 1 else 0.005

    def err_fn(parameters):
        return calculate_error_mj_2d(parameters, time, vel, dt, penalty)[0]

    def fit_fn(parameters):
        return calculate_error_mj_2d(parameters, time, vel, dt, penalty)[1]

    return _decompose_generic(time, vel, num_submovements, lb0, ub0, 4, err_fn, fit_fn, optimizer_constraints)


def decompose_3d(time, vel, num_submovements=None, x_rng=(-5.0, 5.0), y_rng=(0.1, 5.0), z_rng=(-5.0, 5.0), criteria=-np.inf, fitting_constraints=None):
    fitting_constraints = {} if fitting_constraints is None else dict(fitting_constraints)
    optimizer_constraints = remove_non_optimizer_constraints(fitting_constraints)
    constraints = resolve_fitting_constraints(optimizer_constraints)
    penalty = resolve_pair_cancellation_penalty(fitting_constraints)

    time = np.asarray(time, dtype=float)
    vel = np.asarray(vel, dtype=float)

    lb0 = [0, constraints['minDuration'], x_rng[0], y_rng[0], z_rng[0]]
    ub0 = [
        max(time[-1] - constraints['minOnsetSpacing'], constraints['minUpperBoundTime']),
        constraints['maxDuration'],
        x_rng[1],
        y_rng[1],
        z_rng[1],
    ]
    dt = float(time[1] - time[0]) if time.size > 1 else 0.005

    def err_fn(parameters):
        return calculate_error_mj_3d(parameters, time, vel, dt, penalty)[0]

    def fit_fn(parameters):
        return calculate_error_mj_3d(parameters, time, vel, dt, penalty)[1]

    return _decompose_generic(time, vel, num_submovements, lb0, ub0, 5, err_fn, fit_fn, optimizer_constraints)
