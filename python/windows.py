"""Adaptive windowed decomposition APIs for 1D/2D/3D."""

from __future__ import annotations

import time as _wall_time
import numpy as np

from decompose import decompose_1d, decompose_2d, decompose_3d
from min_jerk import minimum_jerk_velocity_1d, minimum_jerk_velocity_2d, minimum_jerk_velocity_3d


def _choose_submovement_index(errors, criteria):
    errors = np.asarray(errors, dtype=float)
    candidates = np.where(errors <= criteria)[0]
    if candidates.size:
        return int(candidates[0])
    candidates = np.where(errors <= 0.05)[0]
    if candidates.size:
        return int(candidates[0])
    candidates = np.where(errors < 0.1)[0]
    if candidates.size:
        return int(candidates[0])
    return int(np.nanargmin(errors))


def _add_component_velocity(dimensions, params_row, t):
    t0 = params_row[0]
    duration = params_row[1]
    amp = params_row[2:]
    if dimensions == 1:
        return minimum_jerk_velocity_1d(t0, duration, amp[0], t)[:, None]
    if dimensions == 2:
        vx, vy, _ = minimum_jerk_velocity_2d(t0, duration, amp[0], amp[1], t)
        return np.column_stack([vx, vy])
    vx, vy, vz, _ = minimum_jerk_velocity_3d(t0, duration, amp[0], amp[1], amp[2], t)
    return np.column_stack([vx, vy, vz])


def _run_windows(
    time,
    vel,
    dimensions,
    submovement_range,
    criteria,
    window_size,
    fitting_constraints,
    decompose_fn,
    pps,
):
    time = np.asarray(time, dtype=float)
    vel = np.asarray(vel, dtype=float)
    if dimensions == 1:
        vel = vel.reshape(-1, 1)

    t0s = []
    durations = []
    amplitudes = []
    endtimes = []
    startwindows = []
    endwindows = []

    best_errors = []
    best_parameters = []
    best_velocity = []

    start = _wall_time.time()

    current_window_start = float(time[0])
    w = 0
    while current_window_start < time[-1]:
        w += 1
        current_window_end = min(current_window_start + window_size, float(time[-1]))
        startwindows.append(current_window_start)
        endwindows.append(current_window_end)

        idx = np.where((time >= current_window_start) & (time <= current_window_end))[0]
        if idx.size == 0:
            break

        t_segment = time[idx]
        v_segment = vel[idx, :].copy()

        for k, endt in enumerate(endtimes):
            if endt > current_window_start:
                v_segment -= _add_component_velocity(dimensions, np.hstack([t0s[k], durations[k], amplitudes[k]]), t_segment)

        local_time = t_segment - t_segment[0]
        errors, parameter_candidates, _ = decompose_fn(
            local_time,
            v_segment if dimensions > 1 else v_segment[:, 0],
            num_submovements=list(submovement_range),
            criteria=criteria,
            fitting_constraints=fitting_constraints,
        )

        chosen = _choose_submovement_index(errors, criteria)
        count = int(list(submovement_range)[chosen])

        params = np.asarray(parameter_candidates[chosen], dtype=float)
        if params.ndim == 0 or np.isnan(params).any():
            accepted = np.empty((0, pps))
            next_window_start = current_window_end
        else:
            params = params.reshape(count, pps)
            params[:, 0] += t_segment[0]
            candidate_endtimes = params[:, 0] + params[:, 1]

            keep = candidate_endtimes <= current_window_end
            deferred = ~keep

            if np.any(deferred):
                next_window_start = float(np.min(params[deferred, 0]))
                if next_window_start <= current_window_start:
                    keep[:] = True
                    next_window_start = current_window_end
                else:
                    endwindows[-1] = next_window_start
            else:
                next_window_start = current_window_end

            accepted = params[keep, :]

        if accepted.shape[0] == 0:
            best_errors.append(np.nan)
            best_parameters.append(np.array([]))
            best_velocity.append(np.zeros((t_segment.size, dimensions if dimensions > 1 else 1)))
        else:
            reconstructed = np.zeros((t_segment.size, dimensions if dimensions > 1 else 1))
            for row in accepted:
                reconstructed += _add_component_velocity(dimensions, row, t_segment)

            denom = float(np.sum(v_segment**2))
            if denom == 0:
                denom = 1.0
            current_error = float(np.sum((reconstructed - v_segment) ** 2) / denom)

            best_errors.append(current_error)
            best_parameters.append(accepted.ravel())
            best_velocity.append(reconstructed)

            for row in accepted:
                t0s.append(float(row[0]))
                durations.append(float(row[1]))
                amplitudes.append(np.asarray(row[2:], dtype=float))
                endtimes.append(float(row[0] + row[1]))

        elapsed = _wall_time.time() - start
        hours = int(elapsed // 3600)
        minutes = int((elapsed - hours * 3600) // 60)
        seconds = int(round(elapsed - hours * 3600 - minutes * 60))
        processed = endwindows[-1] - time[0]
        total = time[-1] - time[0]
        percent = 0.0 if total <= 0 else processed / total * 100
        print(
            f'Finished window {w}, time since start: {hours} hours, {minutes} minutes, {seconds} seconds, '
            f'processed {percent:.1f}% ({processed:.3f} seconds from {total:.3f} seconds)'
        )

        if next_window_start >= time[-1] or processed >= total - 0.01:
            break
        current_window_start = next_window_start

    t0s_arr = np.asarray(t0s, dtype=float)
    d_arr = np.asarray(durations, dtype=float)
    end_arr = np.asarray(endtimes, dtype=float)

    if amplitudes:
        amp_arr = np.vstack(amplitudes)
    else:
        amp_arr = np.empty((0, dimensions))

    submovements_velocity = np.zeros((time.size, len(t0s), dimensions if dimensions > 1 else 1))
    for k in range(len(t0s)):
        row = np.hstack([t0s_arr[k], d_arr[k], amp_arr[k]])
        submovements_velocity[:, k, :] = _add_component_velocity(dimensions, row, time)

    reconstructed = np.sum(submovements_velocity, axis=1)
    if dimensions == 1:
        reconstructed = reconstructed[:, 0]
        submovements_velocity = submovements_velocity[:, :, 0]

    decomposition = {
        't0s': t0s_arr,
        'Ds': d_arr,
        'endtimes': end_arr,
        'time': time,
        'vel': vel[:, 0] if dimensions == 1 else vel,
        'startwindows': np.asarray(startwindows),
        'endwindows': np.asarray(endwindows),
        'submovementsVelocity': submovements_velocity,
        'reconstructedVelocity': reconstructed,
    }

    if dimensions == 1:
        decomposition['As'] = amp_arr[:, 0] if amp_arr.size else np.array([])
        decomposition['parameters'] = np.column_stack([t0s_arr, d_arr, decomposition['As']]).ravel() if t0s_arr.size else np.array([])
    elif dimensions == 2:
        decomposition['Axs'] = amp_arr[:, 0] if amp_arr.size else np.array([])
        decomposition['Ays'] = amp_arr[:, 1] if amp_arr.size else np.array([])
    else:
        decomposition['Axs'] = amp_arr[:, 0] if amp_arr.size else np.array([])
        decomposition['Ays'] = amp_arr[:, 1] if amp_arr.size else np.array([])
        decomposition['Azs'] = amp_arr[:, 2] if amp_arr.size else np.array([])

    return best_errors, best_parameters, best_velocity, decomposition


def decompose_1d_windows(time, vel, submovement_range=(1, 2, 3, 4), a_rng=(-5.0, 5.0), criteria=-np.inf, window_size=3.0, fitting_constraints=None):
    return _run_windows(
        time,
        vel,
        1,
        submovement_range,
        criteria,
        window_size,
        fitting_constraints,
        lambda t, v, num_submovements, criteria, fitting_constraints: decompose_1d(
            t,
            v,
            num_submovements=num_submovements,
            a_rng=a_rng,
            criteria=criteria,
            fitting_constraints=fitting_constraints,
        ),
        3,
    )


def decompose_2d_windows(time, vel, submovement_range=(1, 2, 3, 4), x_rng=(-5.0, 5.0), y_rng=(0.1, 5.0), criteria=-np.inf, window_size=3.0, fitting_constraints=None):
    return _run_windows(
        time,
        vel,
        2,
        submovement_range,
        criteria,
        window_size,
        fitting_constraints,
        lambda t, v, num_submovements, criteria, fitting_constraints: decompose_2d(
            t,
            v,
            num_submovements=num_submovements,
            x_rng=x_rng,
            y_rng=y_rng,
            criteria=criteria,
            fitting_constraints=fitting_constraints,
        ),
        4,
    )


def decompose_3d_windows(
    time,
    vel,
    submovement_range=(1, 2, 3, 4),
    x_rng=(-5.0, 5.0),
    y_rng=(0.1, 5.0),
    z_rng=(-5.0, 5.0),
    criteria=-np.inf,
    window_size=3.0,
    fitting_constraints=None,
):
    return _run_windows(
        time,
        vel,
        3,
        submovement_range,
        criteria,
        window_size,
        fitting_constraints,
        lambda t, v, num_submovements, criteria, fitting_constraints: decompose_3d(
            t,
            v,
            num_submovements=num_submovements,
            x_rng=x_rng,
            y_rng=y_rng,
            z_rng=z_rng,
            criteria=criteria,
            fitting_constraints=fitting_constraints,
        ),
        5,
    )
