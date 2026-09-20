"""Metrics for decomposed submovements."""

from __future__ import annotations

import numpy as np


def calculate_overlap(t0, duration):
    """Compute overlap proportion between all submovement pairs.

    Returns
    -------
    overlap_matrix : ndarray
        Pairwise overlap durations normalized by min duration of each pair.
    mean_overlap : float
        Mean of upper-triangular pairwise overlaps (excluding diagonal).
    """
    t0 = np.asarray(t0, dtype=float)
    duration = np.asarray(duration, dtype=float)
    n = t0.size
    if n == 0:
        return np.zeros((0, 0)), 0.0
    if n == 1:
        return np.zeros((1, 1)), 0.0

    starts = t0[:, None]
    ends = (t0 + duration)[:, None]

    inter_start = np.maximum(starts, starts.T)
    inter_end = np.minimum(ends, ends.T)
    overlap = np.maximum(0.0, inter_end - inter_start)

    min_duration = np.minimum(duration[:, None], duration[None, :])
    with np.errstate(divide='ignore', invalid='ignore'):
        overlap_norm = np.where(min_duration > 0, overlap / min_duration, 0.0)

    np.fill_diagonal(overlap_norm, 0.0)

    triu = np.triu_indices(n, k=1)
    mean_overlap = float(np.mean(overlap_norm[triu])) if triu[0].size else 0.0
    return overlap_norm, mean_overlap


def calculate_relative_onset(t0):
    """Relative onset of each submovement from previous one."""
    t0 = np.asarray(t0, dtype=float)
    if t0.size == 0:
        return np.array([], dtype=float)
    if t0.size == 1:
        return np.array([0.0], dtype=float)
    rel = np.zeros_like(t0)
    rel[1:] = np.diff(t0)
    return rel
