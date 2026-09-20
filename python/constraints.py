"""Shared fitting constraint handling."""

from __future__ import annotations


def resolve_fitting_constraints(overrides: dict | None = None) -> dict:
    defaults = {
        'minOnsetSpacing': 0.167,
        'minDuration': 0.167,
        'maxDuration': 1.0,
        'minUpperBoundTime': 0.1,
        'numRestarts': 20,
        'maxFunEvals': 10**13,
        'maxIter': 5000,
    }

    if overrides is None:
        return defaults.copy()
    if not isinstance(overrides, dict):
        raise ValueError('fittingConstraints must be a dict')

    constraints = defaults.copy()
    for key, value in overrides.items():
        if key not in defaults:
            raise ValueError(f'Unknown fitting constraint: {key}')
        constraints[key] = value

    if constraints['minOnsetSpacing'] <= 0:
        raise ValueError('minOnsetSpacing must be > 0')
    if constraints['minDuration'] <= 0:
        raise ValueError('minDuration must be > 0')
    if constraints['maxDuration'] <= 0:
        raise ValueError('maxDuration must be > 0')
    if constraints['maxDuration'] < constraints['minDuration']:
        raise ValueError('maxDuration must be >= minDuration')
    if constraints['minUpperBoundTime'] <= 0:
        raise ValueError('minUpperBoundTime must be > 0')
    if int(constraints['numRestarts']) < 1:
        raise ValueError('numRestarts must be >= 1')
    if constraints['maxFunEvals'] < 1:
        raise ValueError('maxFunEvals must be >= 1')
    if int(constraints['maxIter']) < 1:
        raise ValueError('maxIter must be >= 1')

    constraints['numRestarts'] = int(constraints['numRestarts'])
    constraints['maxIter'] = int(constraints['maxIter'])

    return constraints


def resolve_pair_cancellation_penalty(overrides: dict | None = None) -> dict:
    settings = {
        'weight': 0.0,
        't0Sigma': 0.05,
        'dSigma': 0.05,
    }
    if overrides is None:
        return settings

    mapping = {
        'pairCancellationPenaltyWeight': 'weight',
        'pairCancellationT0Sigma': 't0Sigma',
        'pairCancellationDSigma': 'dSigma',
    }
    for source_key, target_key in mapping.items():
        if source_key in overrides:
            settings[target_key] = float(overrides[source_key])

    if settings['weight'] < 0:
        raise ValueError('pairCancellationPenaltyWeight must be >= 0')
    if settings['t0Sigma'] <= 0:
        raise ValueError('pairCancellationT0Sigma must be > 0')
    if settings['dSigma'] <= 0:
        raise ValueError('pairCancellationDSigma must be > 0')

    return settings


def remove_non_optimizer_constraints(overrides: dict | None = None) -> dict:
    if overrides is None:
        return {}
    cleaned = dict(overrides)
    cleaned.pop('pairCancellationPenaltyWeight', None)
    cleaned.pop('pairCancellationT0Sigma', None)
    cleaned.pop('pairCancellationDSigma', None)
    return cleaned
