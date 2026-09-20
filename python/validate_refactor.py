"""Sanity checks for Python decomposition refactor."""

from __future__ import annotations

import numpy as np

from constraints import resolve_fitting_constraints
from decompose import decompose_1d, decompose_2d
from metrics import calculate_overlap, calculate_relative_onset
from min_jerk import minimum_jerk_velocity_1d, minimum_jerk_velocity_2d


def _assert(condition, message):
    if not condition:
        raise AssertionError(message)


def test_min_jerk_shapes():
    t = np.linspace(0.0, 1.0, 50)
    v1 = minimum_jerk_velocity_1d(t, 0.1, 0.4, 0.8)
    _assert(v1.shape == t.shape, '1D velocity shape mismatch')

    v2x, v2y, v2t = minimum_jerk_velocity_2d(t, 0.1, 0.4, 1.0, -0.5)
    _assert(v2x.shape == t.shape and v2y.shape == t.shape and v2t.shape == t.shape, '2D velocity shape mismatch')


def test_constraints_resolution():
    c = resolve_fitting_constraints(None)
    _assert('numRestarts' in c and c['numRestarts'] >= 1, 'Invalid default constraints')


def test_decompose_1d_simple():
    t = np.linspace(0.0, 1.0, 200)
    v = minimum_jerk_velocity_1d(t, 0.2, 0.4, 1.2)
    _, p, _ = decompose_1d(t, v, num_submovements=1)
    p = np.asarray(p).reshape(-1, 3)
    _assert(p.ndim == 2 and p.shape[1] == 3, '1D parameters shape invalid')
    _assert(p.shape[0] >= 1, '1D decomposition found no submovements')


def test_decompose_2d_simple():
    t = np.linspace(0.0, 1.0, 200)
    vx, vy, _ = minimum_jerk_velocity_2d(t, 0.2, 0.4, 1.0, 0.3)
    v = np.column_stack((vx, vy))
    _, p, _ = decompose_2d(t, v, num_submovements=1)
    p = np.asarray(p).reshape(-1, 4)
    _assert(p.ndim == 2 and p.shape[1] == 4, '2D parameters shape invalid')
    _assert(p.shape[0] >= 1, '2D decomposition found no submovements')


def test_metrics():
    t0 = np.array([0.1, 0.2, 0.5])
    d = np.array([0.4, 0.3, 0.2])
    overlap, mean_overlap = calculate_overlap(t0, d)
    _assert(overlap.shape == (3, 3), 'Overlap matrix shape invalid')
    _assert(0.0 <= mean_overlap <= 1.0, 'Mean overlap out of range')

    rel = calculate_relative_onset(t0)
    _assert(rel.shape == t0.shape, 'Relative onset shape invalid')
    _assert(np.isclose(rel[0], 0.0), 'First relative onset should be zero')


def main():
    tests = [
        test_min_jerk_shapes,
        test_constraints_resolution,
        test_decompose_1d_simple,
        test_decompose_2d_simple,
        test_metrics,
    ]
    for test in tests:
        test()
        print(f'[OK] {test.__name__}')
    print('All refactor validation tests passed.')


if __name__ == '__main__':
    main()
