"""APF Utils — consolidated utility module (v6.3b).

Merges: _constants.py, _result.py, _linalg.py, _dag.py, _helpers.py
All symbols previously imported from those modules are available here.

Physics files should import:
    from apf.apf_utils import (check, CheckFailure, _result, _zeros, ...)
"""

# ══════════════════════════════════════════════════════════════════════
# SECTION 1: Physical Constants  (was _constants.py)
# ══════════════════════════════════════════════════════════════════════

# Planck 2018 cosmological parameters (arXiv:1807.06209)
PLANCK = {
    'Omega_Lambda': (0.6889, 0.0056),
    'Omega_m':      (0.3111, 0.0056),
    'Omega_b':      (0.0490, 0.0003),
    'Omega_DM':     (0.2621, 0.0056),
    'f_b':          (0.1573, 0.0014),
    'n_s':          (0.9649, 0.0042),
    'r_upper':      0.06,
    'H0':           (67.36, 0.54),
    'sigma8':       (0.8111, 0.0060),
}

# Particle Data Group 2024
PDG = {
    'eta_B':        (6.12e-10, 0.04e-10),
    'alpha_em':     (1/137.035999084, 0),
    'alpha_s':      (0.1181, 0.0011),
    'sin2_theta_W': (0.23122, 0.00003),
    'm_t_pole':     (173.1, 0.9),
    'm_H':          (125.09, 0.24),
    'm_W':          (80.377, 0.012),
    'm_Z':          (91.1876, 0.0021),
}

# Big Bang Nucleosynthesis
BBN = {
    'Y_p':          (0.2449, 0.0040),
    'D_over_H':     (2.547e-5, 0.025e-5),
    'He3_over_H':   (1.1e-5, 0.2e-5),
    'Li7_over_H':   (1.6e-10, 0.3e-10),
}

# Observational bounds (misc)
OBS = {
    'sum_m_nu_upper': 0.120,
    'tau_proton_lower': 1.6e34,
}

# Physical constants
PHYSICAL = {
    'G_N':    6.67430e-11,
    'hbar':   1.054571817e-34,
    'c':      2.99792458e8,
    'k_B':    1.380649e-23,
    'M_Pl':   1.220890e19,
    'v_EW':   246.22,
}


# ══════════════════════════════════════════════════════════════════════
# SECTION 2: Result Builder  (was _result.py)
# ══════════════════════════════════════════════════════════════════════

def result(*, name, tier, epistemic, summary,
           key_result='', dependencies=None, cross_refs=None,
           artifacts=None, **extra):
    """Build a standard theorem result dict.

    Parameters
    ----------
    name : str          Human-readable theorem name
    tier : int          Derivation tier (-1 for axioms, 0+ for theorems)
    epistemic : str     'AXIOM', 'POSTULATE', 'P', 'P_structural', etc.
    summary : str       Full derivation summary
    key_result : str    One-line headline result
    dependencies : list Theorem names this depends on
    cross_refs : list   Related theorem names (not dependencies)
    artifacts : dict    Computed values to store
    **extra             Any additional fields
    """
    out = {
        'passed': True,
        'status': 'PASS',
        'name': name,
        'tier': tier,
        'epistemic': epistemic,
        'summary': summary,
        'key_result': key_result,
        'dependencies': dependencies or [],
        'cross_refs': cross_refs or [],
        'artifacts': artifacts or {},
    }
    out.update(extra)
    return out

# Underscore alias for backward compatibility with all physics modules
_result = result


# ══════════════════════════════════════════════════════════════════════
# SECTION 3: Linear Algebra  (was _linalg.py)
# ══════════════════════════════════════════════════════════════════════

import math as _math


def zeros(n, m=None):
    """n×m zero matrix."""
    if m is None:
        m = n
    return [[0.0] * m for _ in range(n)]


def eye(n):
    """n×n identity."""
    return [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]


def diag(vals):
    """Diagonal matrix from list of values."""
    n = len(vals)
    return [[vals[i] if i == j else 0.0 for j in range(n)] for i in range(n)]


def mat(rows):
    """Copy nested list into matrix."""
    return [list(r) for r in rows]


def mm(A, B):
    """Matrix multiply A @ B."""
    n, p, m = len(A), len(B), len(B[0])
    return [[sum(A[i][k] * B[k][j] for k in range(p))
             for j in range(m)] for i in range(n)]


def madd(A, B):
    """Element-wise A + B."""
    return [[A[i][j] + B[i][j] for j in range(len(A[0]))]
            for i in range(len(A))]


def msub(A, B):
    """Element-wise A - B."""
    return [[A[i][j] - B[i][j] for j in range(len(A[0]))]
            for i in range(len(A))]


def mscale(s, A):
    """Scalar × matrix."""
    return [[s * A[i][j] for j in range(len(A[0]))]
            for i in range(len(A))]


def dag(A):
    """Conjugate transpose (dagger)."""
    n, m = len(A), len(A[0])
    return [[(A[j][i].conjugate() if isinstance(A[j][i], complex)
              else A[j][i]) for j in range(n)] for i in range(m)]


def tr(A):
    """Trace."""
    return sum(A[i][i] for i in range(len(A)))


def det(A):
    """Determinant (up to 4×4, Leibniz expansion for larger)."""
    n = len(A)
    if n == 1:
        return A[0][0]
    if n == 2:
        return A[0][0] * A[1][1] - A[0][1] * A[1][0]
    if n == 3:
        return (A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1])
                - A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0])
                + A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]))
    total = 0
    for j in range(n):
        minor = [[A[r][c] for c in range(n) if c != j]
                 for r in range(1, n)]
        total += ((-1) ** j) * A[0][j] * det(minor)
    return total


def kron(A, B):
    """Kronecker product A ⊗ B."""
    na, ma = len(A), len(A[0])
    nb, mb = len(B), len(B[0])
    return [[A[i // nb][j // mb] * B[i % nb][j % mb]
             for j in range(ma * mb)] for i in range(na * nb)]


def outer(u, v):
    """Outer product |u><v|."""
    return [[u[i] * (v[j].conjugate() if isinstance(v[j], complex) else v[j])
             for j in range(len(v))] for i in range(len(u))]


def partial_trace_B(rho, da, db):
    """Partial trace over subsystem B."""
    return [[sum(rho[i * db + k][j * db + k] for k in range(db))
             for j in range(da)] for i in range(da)]


def vn_entropy(rho):
    """Von Neumann entropy S = -Tr(ρ log ρ)."""
    n = len(rho)
    M = [[rho[i][j].real if isinstance(rho[i][j], complex)
          else float(rho[i][j]) for j in range(n)] for i in range(n)]
    for _ in range(300):
        mx, p, q = 0.0, 0, 1
        for i in range(n):
            for j in range(i + 1, n):
                if abs(M[i][j]) > mx:
                    mx = abs(M[i][j]); p, q = i, j
        if mx < 1e-14:
            break
        if abs(M[p][p] - M[q][q]) < 1e-15:
            theta = _math.pi / 4
        else:
            theta = 0.5 * _math.atan2(2 * M[p][q], M[p][p] - M[q][q])
        c, s = _math.cos(theta), _math.sin(theta)
        Mc = [row[:] for row in M]
        for i in range(n):
            Mc[i][p] = c * M[i][p] + s * M[i][q]
            Mc[i][q] = -s * M[i][p] + c * M[i][q]
        Mr = [row[:] for row in Mc]
        for j in range(n):
            Mr[p][j] = c * Mc[p][j] + s * Mc[q][j]
            Mr[q][j] = -s * Mc[p][j] + c * Mc[q][j]
        M = Mr
    evals = [max(M[i][i], 1e-300) for i in range(n)]
    return -sum(e * _math.log(e) for e in evals if e > 1e-15)


def aclose(a, b, tol=1e-9):
    """Check if a ≈ b elementwise."""
    if isinstance(a, (list, tuple)):
        return all(aclose(ai, bi, tol) for ai, bi in zip(a, b))
    return abs(a - b) < tol


# Underscore aliases expected by physics modules
_zeros = zeros
_eye   = eye
_diag  = diag
_mat   = mat
_mm    = mm
_madd  = madd
_msub  = msub
_mscale = mscale
_dag   = dag
_tr    = tr
_det   = det
_kron  = kron
_outer = outer
_partial_trace_B = partial_trace_B
_vn_entropy = vn_entropy
_aclose = aclose


# ══════════════════════════════════════════════════════════════════════
# SECTION 4: Derivation DAG Cache  (was _dag.py)
# ══════════════════════════════════════════════════════════════════════

from fractions import Fraction as _Fraction

__all__ = [
    # Constants
    'PLANCK', 'PDG', 'BBN', 'OBS', 'PHYSICAL',
    # Result builder
    'result', '_result',
    # Linalg
    'zeros', 'eye', 'diag', 'mat', 'mm', 'madd', 'msub', 'mscale',
    'dag', 'tr', 'det', 'kron', 'outer', 'partial_trace_B', 'vn_entropy', 'aclose',
    '_zeros', '_eye', '_diag', '_mat', '_mm', '_madd', '_msub', '_mscale',
    '_dag', '_tr', '_det', '_kron', '_outer', '_partial_trace_B', '_vn_entropy', '_aclose',
    # DAG
    'dag_put', 'dag_get', 'dag_reset', 'dag_dump', 'dag_has', 'dag_verify_chain',
    'ChainInconsistency',
    # Helpers
    'check', 'CheckFailure',
    '_eigvalsh', '_eigh', '_eigh_3x3',
    '_fnorm', '_mv', '_zvec', '_vdot', '_vkron', '_vscale', '_vadd',
]


class ChainInconsistency(Exception):
    """Raised when a DAG-linked value doesn't match expectations."""
    pass


class _DAGEntry:
    __slots__ = ('key', 'value', 'source', 'derivation', 'consumers')

    def __init__(self, key, value, source, derivation=""):
        self.key = key
        self.value = value
        self.source = source
        self.derivation = derivation
        self.consumers = []

    def __repr__(self):
        return (f"DAGEntry({self.key}={self.value}, "
                f"from={self.source}, used_by={self.consumers})")


def _values_equal(a, b):
    if a is b:
        return True
    if a == b:
        return True
    try:
        fa, fb = float(a), float(b)
        if abs(fa - fb) < 1e-12 * (1 + abs(fa)):
            return True
    except (TypeError, ValueError):
        pass
    return str(a) == str(b)


class _DAGCache:
    def __init__(self):
        self._store = {}
        self._enabled = True

    def reset(self):
        self._store.clear()

    def put(self, key, value, source="?", derivation=""):
        if key in self._store:
            existing = self._store[key]
            if existing.source == source and _values_equal(existing.value, value):
                return
            if not _values_equal(existing.value, value):
                import warnings
                warnings.warn(
                    f"DAG: key '{key}' overwritten by {source} "
                    f"(was {existing.value} from {existing.source}, now {value})",
                    RuntimeWarning, stacklevel=3,
                )
        self._store[key] = _DAGEntry(key, value, source, derivation)

    def get(self, key, default=None, consumer="?", expected_source=None, verify=True):
        if key in self._store:
            entry = self._store[key]
            entry.consumers.append(consumer)
            if verify and default is not None:
                if not _values_equal(entry.value, default):
                    raise ChainInconsistency(
                        f"Chain inconsistency in {consumer}: "
                        f"DAG['{key}'] = {entry.value} (from {entry.source}) "
                        f"but downstream expected {default}."
                    )
            if expected_source and entry.source != expected_source:
                import warnings
                warnings.warn(
                    f"DAG: {consumer} expected '{key}' from "
                    f"{expected_source}, got it from {entry.source}",
                    RuntimeWarning, stacklevel=3,
                )
            return entry.value
        return default

    def has(self, key):
        return key in self._store

    def dump(self):
        return {
            key: {
                'value': entry.value,
                'source': entry.source,
                'derivation': entry.derivation,
                'consumers': entry.consumers,
            }
            for key, entry in self._store.items()
        }

    def verify_chain(self, chain_spec):
        results = []
        for upstream_key, downstream, expected in chain_spec:
            entry = self._store.get(upstream_key)
            if entry is None:
                results.append({'link': f'{upstream_key} → {downstream}',
                                'ok': False, 'detail': f'Key {upstream_key} not in DAG'})
            elif not _values_equal(entry.value, expected):
                results.append({'link': f'{upstream_key} → {downstream}',
                                'ok': False,
                                'detail': f'Expected {expected}, got {entry.value} (from {entry.source})'})
            else:
                consumed = downstream in entry.consumers
                results.append({'link': f'{upstream_key} → {downstream}',
                                'ok': True, 'consumed': consumed,
                                'detail': f'{entry.value} from {entry.source}'})
        return results


_cache = _DAGCache()


def dag_put(key, value, source="?", derivation=""):
    """Store a derived value in the global DAG cache."""
    _cache.put(key, value, source, derivation)


def dag_get(key, default=None, consumer="?", expected_source=None, verify=True):
    """Retrieve a derived value from the global DAG cache."""
    return _cache.get(key, default, consumer, expected_source, verify)


def dag_reset():
    """Clear all cached values."""
    _cache.reset()


def dag_dump():
    """Return full cache contents."""
    return _cache.dump()


def dag_has(key):
    """Check if a key exists in the cache."""
    return _cache.has(key)


def dag_verify_chain(chain_spec):
    """Verify a derivation chain."""
    return _cache.verify_chain(chain_spec)


# ══════════════════════════════════════════════════════════════════════
# SECTION 5: Helpers & Eigensolvers  (was _helpers.py)
# ══════════════════════════════════════════════════════════════════════


class CheckFailure(Exception):
    """Raised when a theorem check fails (never stripped by python -O)."""
    pass


def check(condition, msg=""):
    """Assert replacement immune to python -O."""
    if not condition:
        raise CheckFailure(msg)


def _eigvalsh(A):
    """Eigenvalues of Hermitian matrix (Jacobi iteration, real projection).
    Backward-compatible with v4.x monolith behavior.
    """
    n = len(A)
    M = [[A[i][j].real if isinstance(A[i][j], complex) else float(A[i][j])
          for j in range(n)] for i in range(n)]
    for _ in range(300):
        p, q, mx = 0, 1, 0.0
        for i in range(n):
            for j in range(i + 1, n):
                if abs(M[i][j]) > mx:
                    mx = abs(M[i][j]); p, q = i, j
        if mx < 1e-14:
            break
        if abs(M[p][p] - M[q][q]) < 1e-15:
            theta = _math.pi / 4
        else:
            theta = 0.5 * _math.atan2(2 * M[p][q], M[p][p] - M[q][q])
        c, s = _math.cos(theta), _math.sin(theta)
        Mc = [row[:] for row in M]
        for i in range(n):
            Mc[i][p] = c * M[i][p] + s * M[i][q]
            Mc[i][q] = -s * M[i][p] + c * M[i][q]
        Mr = [row[:] for row in Mc]
        for j in range(n):
            Mr[p][j] = c * Mc[p][j] + s * Mc[q][j]
            Mr[q][j] = -s * Mc[p][j] + c * Mc[q][j]
        M = Mr
    return sorted(M[i][i] for i in range(n))


def _zvec(n):
    """Zero vector of length n."""
    return [complex(0)] * n


def _mv(A, v):
    """Matrix-vector multiply."""
    return [sum(A[i][k] * v[k] for k in range(len(v))) for i in range(len(A))]


def _fnorm(A):
    """Frobenius norm."""
    return _math.sqrt(sum(abs(A[i][j]) ** 2
                          for i in range(len(A))
                          for j in range(len(A[0]))))


def _vdot(u, v):
    """Complex inner product <u|v>."""
    return sum(a.conjugate() * b for a, b in zip(u, v))


def _vkron(u, v):
    """Kronecker product of vectors."""
    return [a * b for a in u for b in v]


def _vscale(c, v):
    """Scalar * vector."""
    c = complex(c)
    return [c * x for x in v]


def _vadd(u, v):
    """Vector addition."""
    return [a + b for a, b in zip(u, v)]


def _eigh_3x3(H):
    """Eigenvalues and eigenvectors for 3×3 Hermitian. Delegates to _eigh."""
    return _eigh(H)


def _eigh(H):
    """Eigenvalues and eigenvectors of Hermitian matrix via two-step Jacobi."""
    n = len(H)
    A = [[complex(H[i][j]) for j in range(n)] for i in range(n)]
    V = [[complex(1 if i == j else 0) for j in range(n)] for i in range(n)]

    for _ in range(800):
        mx, p, q = 0.0, 0, 1
        for i in range(n):
            for j in range(i + 1, n):
                if abs(A[i][j]) > mx:
                    mx = abs(A[i][j]); p, q = i, j
        if mx < 1e-13:
            break
        apq = A[p][q]
        r = abs(apq)
        if r < 1e-15:
            continue
        eia = apq / r
        eia_c = eia.conjugate()
        for i in range(n):
            A[i][q] *= eia_c
        for j in range(n):
            A[q][j] *= eia
        A[p][q] = complex(r)
        A[q][p] = complex(r)
        for i in range(n):
            V[i][q] *= eia_c
        app, aqq = A[p][p].real, A[q][q].real
        if abs(app - aqq) < 1e-15:
            theta = _math.pi / 4
        else:
            theta = 0.5 * _math.atan2(2 * r, app - aqq)
        c, s = _math.cos(theta), _math.sin(theta)
        for i in range(n):
            aip, aiq = A[i][p], A[i][q]
            A[i][p] = c * aip + s * aiq
            A[i][q] = -s * aip + c * aiq
        for j in range(n):
            apj, aqj = A[p][j], A[q][j]
            A[p][j] = c * apj + s * aqj
            A[q][j] = -s * apj + c * aqj
        A[p][p] = complex(A[p][p].real)
        A[q][q] = complex(A[q][q].real)
        A[p][q] = complex(0)
        A[q][p] = complex(0)
        for i in range(n):
            vip, viq = V[i][p], V[i][q]
            V[i][p] = c * vip + s * viq
            V[i][q] = -s * vip + c * viq
    evals = [A[i][i].real for i in range(n)]
    idx = sorted(range(n), key=lambda i: evals[i])
    return ([evals[i] for i in idx],
            [[V[r][idx[c]] for c in range(n)] for r in range(n)])
