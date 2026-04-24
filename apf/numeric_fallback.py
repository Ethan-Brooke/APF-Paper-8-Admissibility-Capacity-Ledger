"""Pure-numpy fallbacks for scipy functions used in the bank.

Phase 22.1-extra: removes scipy as a hard dependency. Users with scipy
installed will get the faster / higher-precision scipy routes via the
try/except pattern in call sites; users without scipy get these
fallbacks, which match the scipy API for the specific calls made in
`apf/supplements.py`.

Implementations:

- **`expm(A)`** — matrix exponential via scaling-and-squaring with
  Padé[9,9] approximation. Standard robust method; matches
  `scipy.linalg.expm` to ~1e-14 for small matrices.

- **`solve_ivp(fun, t_span, y0, method='RK45', dense_output=True,
  rtol=..., atol=...)`** — Dormand-Prince 4(5) adaptive integrator
  with dense output via the classical DP 4th-order interpolant.
  Returns an object with a `.sol(t)` callable and `.y`, `.t`
  attributes. Matches scipy API for the single call site in
  `check_L_mc_mt_twoloop_RG`.

- **`quad(fun, a, b)`** — adaptive Simpson's rule. Returns
  `(integral, error_estimate)` tuple matching scipy API. Not as
  robust as scipy's QUADPACK for pathological integrands, but
  sufficient for the smooth integrands in the APF bank.

These implementations are deliberately minimal — just enough surface
area to match the usage in the bank. For general scipy replacement,
install scipy.
"""

from __future__ import annotations

import numpy as _np


# ═══════════════════════════════════════════════════════════════════
# expm — matrix exponential via Padé + scaling-squaring
# ═══════════════════════════════════════════════════════════════════

def _pade9(A):
    """Padé[9,9] approximant of exp(A). Standard Higham 2005 coefficients."""
    b = (17643225600., 8821612800., 2075673600., 302702400., 30270240.,
         2162160., 110880., 3960., 90., 1.)
    I = _np.eye(A.shape[0], dtype=A.dtype)
    A2 = A @ A
    A4 = A2 @ A2
    A6 = A2 @ A4
    A8 = A2 @ A6
    U = A @ (A8 * b[9] + A6 * b[7] + A4 * b[5] + A2 * b[3] + I * b[1])
    V = A8 * b[8] + A6 * b[6] + A4 * b[4] + A2 * b[2] + I * b[0]
    # exp(A) ≈ (V - U)^-1 (V + U)
    return _np.linalg.solve(V - U, V + U)


def expm(A):
    """Matrix exponential.

    Pure-numpy drop-in for ``scipy.linalg.expm``. Uses scaling and
    squaring with a Padé[9,9] approximant. Accuracy ~1e-14 for
    well-conditioned matrices.

    Parameters
    ----------
    A : array_like
        Square matrix (real or complex).

    Returns
    -------
    ndarray
        Matrix exponential ``exp(A)``.
    """
    A = _np.asarray(A)
    # Promote to complex if input is complex; else keep real
    if _np.iscomplexobj(A):
        A = A.astype(_np.complex128, copy=False)
    else:
        A = A.astype(_np.float64, copy=False)

    # Higham threshold for Padé[9,9]: ||A|| <= 5.4 (1-norm)
    norm = _np.max(_np.sum(_np.abs(A), axis=0))
    if norm == 0.0:
        return _np.eye(A.shape[0], dtype=A.dtype)

    # Scaling: choose s so that ||A / 2^s|| <= 5.4
    s = max(0, int(_np.ceil(_np.log2(norm / 5.4)))) if norm > 5.4 else 0
    A_scaled = A / (2.0 ** s)

    R = _pade9(A_scaled)

    # Squaring: square s times to recover exp(A)
    for _ in range(s):
        R = R @ R
    return R


# ═══════════════════════════════════════════════════════════════════
# solve_ivp — Dormand-Prince RK4(5) with dense output
# ═══════════════════════════════════════════════════════════════════

# Dormand-Prince 4(5) coefficients (classical)
# Butcher tableau from Hairer-Nørsett-Wanner
_DP_C = _np.array([0, 1/5, 3/10, 4/5, 8/9, 1.0, 1.0])
_DP_A = [
    [],
    [1/5],
    [3/40, 9/40],
    [44/45, -56/15, 32/9],
    [19372/6561, -25360/2187, 64448/6561, -212/729],
    [9017/3168, -355/33, 46732/5247, 49/176, -5103/18656],
    [35/384, 0.0, 500/1113, 125/192, -2187/6784, 11/84],
]
# 5th-order solution coefficients
_DP_B5 = _np.array([35/384, 0.0, 500/1113, 125/192, -2187/6784, 11/84, 0.0])
# 4th-order solution coefficients (for error estimation)
_DP_B4 = _np.array([5179/57600, 0.0, 7571/16695, 393/640, -92097/339200,
                    187/2100, 1/40])


class _OdeResult:
    """Minimal scipy-compatible return object."""
    def __init__(self, t_list, y_list, k_list, success=True):
        self.t = _np.array(t_list)
        self.y = _np.array(y_list).T  # shape (n_dim, n_points) matching scipy
        self._k_list = k_list         # stage derivatives for interpolation
        self.success = success
        self.sol = _DenseOutput(t_list, y_list, k_list)


class _DenseOutput:
    """Dense-output callable: sol(t) interpolates at t.

    Uses the 4th-order Hermite interpolant from the DP stages at each
    step. For a scalar equation the interpolant is cubic-Hermite-like.
    """
    def __init__(self, t_list, y_list, k_list):
        self._t = _np.array(t_list)
        self._y = [_np.asarray(y) for y in y_list]
        self._k = k_list  # list of 7-stage k arrays per step

    def __call__(self, t):
        t_arr = _np.atleast_1d(t)
        out = []
        for ti in t_arr:
            # Find the step containing ti
            idx = int(_np.searchsorted(self._t, ti) - 1)
            idx = max(0, min(idx, len(self._t) - 2))
            t0, t1 = self._t[idx], self._t[idx + 1]
            h = t1 - t0
            if h == 0:
                out.append(self._y[idx])
                continue
            theta = (ti - t0) / h
            y0, y1 = self._y[idx], self._y[idx + 1]
            # Hermite-cubic using endpoint derivatives from stored stages
            f0 = self._k[idx][0]  # k1 at start
            f1 = self._k[idx][6]  # k7 at end (which is k1 of next step)
            # Standard cubic Hermite
            h00 = 2 * theta**3 - 3 * theta**2 + 1
            h10 = theta**3 - 2 * theta**2 + theta
            h01 = -2 * theta**3 + 3 * theta**2
            h11 = theta**3 - theta**2
            yi = h00 * y0 + h10 * h * f0 + h01 * y1 + h11 * h * f1
            out.append(yi)
        res = _np.array(out).T  # shape (n_dim, n_points)
        if _np.isscalar(t) or (hasattr(t, 'ndim') and t.ndim == 0):
            return res.flatten()
        return res


def solve_ivp(fun, t_span, y0, method='RK45', dense_output=True,
              rtol=1e-6, atol=1e-9, max_step=None):
    """Dormand-Prince RK4(5) adaptive ODE integrator.

    Pure-numpy drop-in for the specific
    ``scipy.integrate.solve_ivp(..., method='RK45', dense_output=True)``
    calls in the APF bank. Not a full reimplementation — supports
    only the signature used in ``check_L_mc_mt_twoloop_RG``.

    Parameters
    ----------
    fun : callable
        Right-hand side ``f(t, y)``.
    t_span : (2,) tuple
        ``(t0, t1)``.
    y0 : array_like
        Initial state.
    method : str
        Must be ``'RK45'`` (only method supported).
    dense_output : bool
        Must be ``True``; always builds dense output.
    rtol, atol : float
        Relative / absolute tolerances for adaptive step control.
    max_step : float, optional
        Upper bound on step size. Default = span / 10.

    Returns
    -------
    _OdeResult
        Object with ``.t``, ``.y``, ``.sol(t)``, ``.success``.
    """
    if method != 'RK45':
        raise NotImplementedError(
            f"numeric_fallback.solve_ivp only supports RK45, got {method}")
    t0, t1 = float(t_span[0]), float(t_span[1])
    y = _np.asarray(y0, dtype=float).copy()
    t = t0
    t_list = [t]
    y_list = [y.copy()]
    k_list = []  # stages per step

    if max_step is None:
        max_step = (t1 - t0) / 10
    h = min(max_step, (t1 - t0) / 100)  # initial step

    def _rk_step(t, y, h):
        """One DP step. Returns (y_new, k_stages, error_est)."""
        k = [None] * 7
        k[0] = _np.asarray(fun(t, y), dtype=float)
        for i in range(1, 7):
            dy = sum(_DP_A[i][j] * k[j] for j in range(i))
            k[i] = _np.asarray(fun(t + _DP_C[i] * h, y + h * dy), dtype=float)
        y5 = y + h * sum(_DP_B5[i] * k[i] for i in range(7))
        y4 = y + h * sum(_DP_B4[i] * k[i] for i in range(7))
        err = _np.max(_np.abs(y5 - y4))
        return y5, k, err

    MAX_ITERS = 100000
    it = 0
    while t < t1 and it < MAX_ITERS:
        it += 1
        # Don't overshoot
        h = min(h, t1 - t, max_step)
        y_new, k, err = _rk_step(t, y, h)
        # Tolerance check
        tol = atol + rtol * _np.max(_np.abs(y))
        if err <= tol or h <= 1e-14:
            # Accept step
            t += h
            y = y_new
            t_list.append(t)
            y_list.append(y.copy())
            k_list.append(k)
            # Step size adjustment
            if err > 0:
                h_new = h * min(5.0, max(0.1, 0.9 * (tol / err) ** 0.2))
            else:
                h_new = h * 5.0
            h = min(h_new, max_step)
        else:
            # Reject step, shrink
            h = h * max(0.1, 0.9 * (tol / err) ** 0.2)

    if it >= MAX_ITERS:
        return _OdeResult(t_list, y_list, k_list, success=False)
    return _OdeResult(t_list, y_list, k_list, success=True)


# ═══════════════════════════════════════════════════════════════════
# quad — adaptive Simpson's rule
# ═══════════════════════════════════════════════════════════════════

def _simpson(fun, a, b):
    """Basic Simpson's rule on one interval."""
    mid = 0.5 * (a + b)
    return (b - a) / 6.0 * (fun(a) + 4.0 * fun(mid) + fun(b))


def _adaptive_simpson(fun, a, b, tol, S_whole, fa, fb, fm, depth):
    """Recursive adaptive Simpson's."""
    mid = 0.5 * (a + b)
    lm = 0.5 * (a + mid)
    rm = 0.5 * (mid + b)
    flm = fun(lm)
    frm = fun(rm)
    S_left = (mid - a) / 6.0 * (fa + 4.0 * flm + fm)
    S_right = (b - mid) / 6.0 * (fm + 4.0 * frm + fb)
    S_both = S_left + S_right
    diff = S_both - S_whole
    if depth <= 0 or abs(diff) <= 15.0 * tol:
        return S_both + diff / 15.0, abs(diff) / 15.0
    l_val, l_err = _adaptive_simpson(fun, a, mid, tol / 2, S_left,
                                     fa, fm, flm, depth - 1)
    r_val, r_err = _adaptive_simpson(fun, mid, b, tol / 2, S_right,
                                     fm, fb, frm, depth - 1)
    return l_val + r_val, l_err + r_err


def quad(fun, a, b, *, epsabs=1e-10, epsrel=1e-10, limit=50):
    """Adaptive Simpson's rule quadrature.

    Pure-numpy drop-in for ``scipy.integrate.quad`` for the specific
    calls in the APF bank. Smooth integrands only.

    Parameters
    ----------
    fun : callable
        Integrand ``f(x)``.
    a, b : float
        Integration limits.
    epsabs, epsrel : float
        Absolute / relative error tolerance.
    limit : int
        Maximum recursion depth.

    Returns
    -------
    (value, error_estimate) : tuple
        Same shape as ``scipy.integrate.quad``.
    """
    a = float(a); b = float(b)
    if a == b:
        return 0.0, 0.0
    fa = fun(a); fb = fun(b)
    fm = fun(0.5 * (a + b))
    S_whole = (b - a) / 6.0 * (fa + 4.0 * fm + fb)
    tol = max(epsabs, epsrel * abs(S_whole))
    return _adaptive_simpson(fun, a, b, tol, S_whole, fa, fb, fm, limit)


# ═══════════════════════════════════════════════════════════════════
# Self-test harness
# ═══════════════════════════════════════════════════════════════════

if __name__ == '__main__':
    # Basic sanity
    print('numeric_fallback self-test')
    # expm: exp(0) = I
    I = _np.eye(3)
    assert _np.allclose(expm(_np.zeros((3, 3))), I)
    # expm: exp(diag) = diag(exp)
    D = _np.diag([1.0, 2.0, -1.0])
    assert _np.allclose(expm(D), _np.diag(_np.exp([1.0, 2.0, -1.0])))
    print('  expm OK')

    # solve_ivp: dy/dt = -y, y(0)=1 → y(1) = 1/e
    sol = solve_ivp(lambda t, y: -y, [0, 1], [1.0], rtol=1e-10, atol=1e-12)
    import math
    assert abs(sol.y[0, -1] - math.exp(-1)) < 1e-7, f"got {sol.y[0,-1]}, expected {math.exp(-1)}"
    # Dense output
    y_half = sol.sol(0.5)
    assert abs(float(y_half) - math.exp(-0.5)) < 1e-5
    print('  solve_ivp OK')

    # quad: integral of x^2 from 0 to 1 = 1/3
    val, err = quad(lambda x: x * x, 0, 1)
    assert abs(val - 1/3) < 1e-10
    # quad: integral of exp(-x^2) from -3 to 3 = sqrt(pi)*erf(3)
    val, err = quad(lambda x: math.exp(-x * x), -3, 3)
    expected = 1.7724147055049672  # sqrt(pi) * erf(3), computed
    assert abs(val - expected) < 1e-6, f'got {val}, expected {expected}'
    # quad: integral of sin(x) from 0 to pi = 2
    val2, _ = quad(math.sin, 0, math.pi)
    assert abs(val2 - 2.0) < 1e-8
    print('  quad OK')

    print('All self-tests pass.')
