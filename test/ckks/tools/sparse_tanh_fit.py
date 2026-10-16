"""Offline bounded-L1 sine fit; needs numpy/scipy, no keys or ciphertexts.

OPENBLAS_NUM_THREADS=1 python3 test/ckks/tools/sparse_tanh_fit.py
The separate C++ sparse_tanh_certificate.cpp checks the published dyadic fit
without relying on the optimizer or its sampled LP objective.
"""
import json
import numpy as np
import scipy
from scipy.optimize import linprog


def closure(roots):
    return sorted(set(roots) | {a+b for a in roots for b in roots if a+b <= 40}
                  | {abs(a-b) for a in roots for b in roots if a != b})


def fit(roots=(7, 9, 10, 15), norm_cap=2.0):
    frequencies = closure(roots)
    x = np.linspace(0, 4, 2049)[1:]
    target = np.tanh(x)
    matrix = np.sin(x[:, None]*np.array(frequencies)[None, :]*np.pi/8)
    count = len(frequencies)
    # Split sine coefficients into nonnegative parts to bound their L1 norm.
    # Scale only the residual variable, retaining well-scaled coefficients.
    scale = 1e-7
    constraints = np.vstack([
        np.c_[matrix, -matrix, -scale*np.ones(len(x))],
        np.c_[-matrix, matrix, -scale*np.ones(len(x))],
        np.r_[np.ones(2*count), 0][None, :]])
    result = linprog(np.r_[np.zeros(2*count), 1.], A_ub=constraints,
                     b_ub=np.r_[target, -target, norm_cap],
                     bounds=[(0, None)]*(2*count+1), method='highs-ds',
                     options={'dual_feasibility_tolerance': 1e-10,
                              'primal_feasibility_tolerance': 1e-10})
    if not result.success:
        raise RuntimeError(result.message)
    coefficients = result.x[:count]-result.x[count:2*count]
    # Complex series weights c_k=-i*a_k/2 at the evaluator's 28-bit scale.
    weights = np.rint(-coefficients*2**27).astype(np.int64)
    rounded = -weights/2**27
    sampled = 0.0
    intervals = 1 << 20
    for first in range(0, intervals+1, 8192):
        points = np.arange(first, min(first+8192, intervals+1))*4/intervals
        values = np.sin(points[:, None]*np.array(frequencies)[None, :]*np.pi/8)@rounded
        sampled = max(sampled, float(np.max(np.abs(values-np.tanh(points)))))
    curvature = 2+float(np.sum(np.abs(rounded)*(np.array(frequencies)*np.pi/8)**2))
    between_grid = curvature*(4/intervals)**2/8
    return {'roots': list(roots), 'frequencies': frequencies,
            'integer_imaginary_weights28': weights.tolist(),
            'scipy_version': scipy.__version__, 'numpy_version': np.__version__,
            'norm_cap': norm_cap, 'rounded_l1': float(np.abs(rounded).sum()),
            'lp_sampled_objective': float(result.fun*scale),
            'rounded_sampled_error': sampled, 'curvature_bound': curvature,
            'between_grid_bound': between_grid, 'floating_allowance': 1e-12,
            'uniform_bound': sampled+between_grid+1e-12,
            'qualification': 'Public approximation only; encrypted precision unproven.'}


if __name__ == '__main__':
    print(json.dumps(fit(), indent=2))
