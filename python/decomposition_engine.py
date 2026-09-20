"""Shared decomposition optimizer engine for 1D/2D/3D wrappers."""

from __future__ import annotations

import numpy as np
from scipy.optimize import minimize

from constraints import resolve_fitting_constraints


def decompose_nd(
    time: np.ndarray,
    num_submovements: int,
    lower_bound_per_sub: np.ndarray,
    upper_bound_per_sub: np.ndarray,
    parameters_per_submovement: int,
    error_function,
    fit_function,
    fitting_constraints: dict | None = None,
):
    constraints = resolve_fitting_constraints(fitting_constraints)

    if time.size == 0:
        return np.nan, np.full(num_submovements * parameters_per_submovement, np.nan), np.nan

    if np.any(lower_bound_per_sub > upper_bound_per_sub):
        raise ValueError('Lower bounds exceed upper bound - infeasible')

    lb = np.empty(num_submovements * parameters_per_submovement, dtype=float)
    ub = np.empty(num_submovements * parameters_per_submovement, dtype=float)

    for i in range(num_submovements):
        current_lb = lower_bound_per_sub.copy()
        current_lb[0] = i * constraints['minOnsetSpacing']
        if current_lb[0] > upper_bound_per_sub[0]:
            return np.nan, np.nan, np.nan

        start = i * parameters_per_submovement
        end = (i + 1) * parameters_per_submovement
        lb[start:end] = current_lb
        ub[start:end] = upper_bound_per_sub

    bounds = list(zip(lb, ub))

    best_error = np.inf
    best_parameters = None

    for _ in range(constraints['numRestarts']):
        init_params = np.empty(num_submovements * parameters_per_submovement, dtype=float)
        for i in range(num_submovements):
            start = i * parameters_per_submovement
            end = (i + 1) * parameters_per_submovement
            init_params[start:end] = lower_bound_per_sub + (upper_bound_per_sub - lower_bound_per_sub) * np.random.rand(parameters_per_submovement)

        try:
            result = minimize(
                lambda p: float(error_function(p)),
                x0=init_params,
                method='L-BFGS-B',
                bounds=bounds,
                options={
                    'maxiter': constraints['maxIter'],
                    'maxfun': int(min(constraints['maxFunEvals'], 10**9)),
                },
            )
            epsilon = float(error_function(result.x))
            if np.isfinite(epsilon) and epsilon < best_error and np.isrealobj(result.x):
                best_error = epsilon
                best_parameters = result.x
        except Exception:
            continue

    if best_parameters is None:
        return np.nan, np.nan, np.nan

    parameter_matrix = best_parameters.reshape(num_submovements, parameters_per_submovement)
    parameter_matrix = parameter_matrix[np.argsort(parameter_matrix[:, 0]), :]
    best_parameters = parameter_matrix.ravel()

    best_fit = fit_function(best_parameters)

    return best_error, best_parameters, best_fit
