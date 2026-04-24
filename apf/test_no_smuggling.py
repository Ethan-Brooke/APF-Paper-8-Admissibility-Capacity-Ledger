"""Anti-smuggling tests — Phase 22 code-level audit discipline.

The concern addressed by this file, stated bluntly: verifying that
{3, 16, 42, 61, 102} behave consistently AFTER those numbers have been
installed into an ACC record or an operator model is not the same as
deriving them. A verification appendix can contain such checks
legitimately, but a *derivation* claim must be defended against the
replacement attack: if you swap the installed number for a different
one, does the structural check fail for a structural reason, or does
it merely flag a hardcoded-equality disagreement?

These tests probe the boundary. A passing test means: the structural
machinery rejects the mutation for reasons upstream of the
replacement, so the original passing check was not begging the
question.

Organising principle: each test names (a) the claim under attack, (b)
the specific mutation applied, (c) the structural reason the mutation
should fail, (d) the observed behaviour. A test that merely asserts
"mutated value != original value" is a tautology and is not welcome
here.

Scope: this file is a pilot for Phase 22, sized for the acc_SM
residual-provenance fix that landed 2026-04-24. Corpus-wide rollout
(mutation tests on every structural theorem, claim_type metadata tags,
smuggling audit table in the supplement) follows in subsequent
sessions.

Status: every test below is expected to PASS. A failing test here is
a red flag that a smuggling hole was introduced or was always there.

Usage
-----
Run standalone::

    python3 -m apf.test_no_smuggling

Or as part of a future automated suite. These tests are deliberately
NOT registered with the bank REGISTRY — they are adversarial, they
mutate state, and they run under verify_all only through the
top-level guard at the bottom of this file.
"""

from fractions import Fraction

from apf.apf_utils import dag_get, dag_put, dag_has, dag_reset
from apf.unification import (
    ACC, acc_SM,
    _CANON_OMEGA_B_NUM, _CANON_OMEGA_C_NUM, _CANON_OMEGA_LAMBDA_NUM,
    _identity_I2_gauge_cosmological,
)


# ──────────────────────────────────────────────────────────────────
# Test 1: acc_SM honours upstream DAG values when L_equip has run
# ──────────────────────────────────────────────────────────────────

def test_acc_SM_reads_residuals_from_dag_when_populated():
    """Claim under attack: acc_SM() hardcodes the residual partition.

    Mutation: write deliberately non-canonical residual sector counts
    into the DAG (n_baryon=4, n_dark=15, n_vacuum=42 — still summing
    to 61) under a made-up source tag, then call acc_SM() and inspect
    the returned record's residuals.

    Structural expectation: the returned ACC record reflects the
    DAG values, not the CANON constants. If acc_SM() is truly
    provenance-driven, swapping the upstream values must swap the
    ACC record.

    Observed behaviour required: PASS iff acc.residuals == the
    mutated DAG values.
    """
    # Set up a clean DAG state with the three upstream keys that
    # acc_SM consumes, PLUS the K/d_eff keys it already consumes.
    dag_reset()
    dag_put('C_total', 61, source='L_count',
            derivation='test: 45+4+12')
    dag_put('d_eff', 102, source='L_self_exclusion',
            derivation='test: (C_total-1)+C_vacuum')
    # Mutated residuals (still sum to 61; no tautology).
    dag_put('n_baryon', 4, source='test_smuggling',
            derivation='mutation: 3 -> 4 to probe provenance')
    dag_put('n_dark', 15, source='test_smuggling',
            derivation='mutation: 16 -> 15 to probe provenance')
    dag_put('n_vacuum', 42, source='test_smuggling',
            derivation='mutation: unchanged (still structural 42)')

    acc = acc_SM()

    assert acc.residuals['Omega_b'] == 4, (
        f"acc_SM ignored DAG mutation. Expected Omega_b=4 (from DAG), "
        f"got {acc.residuals['Omega_b']}. acc_SM is still hardcoding.")
    assert acc.residuals['Omega_c'] == 15, (
        f"acc_SM ignored DAG mutation. Expected Omega_c=15 (from DAG), "
        f"got {acc.residuals['Omega_c']}.")
    assert acc.residuals['Omega_Lambda'] == 42, (
        f"acc_SM ignored DAG mutation. Expected Omega_Lambda=42 (from DAG), "
        f"got {acc.residuals['Omega_Lambda']}.")

    dag_reset()
    return {'passed': True, 'claim': 'acc_SM is provenance-driven'}


# ──────────────────────────────────────────────────────────────────
# Test 2: I2 consistency passes on DAG-sourced residuals (arithmetic
# identity is what it is — this is the smallest possible witness
# that the fix doesn't regress I2)
# ──────────────────────────────────────────────────────────────────

def test_I2_passes_on_canonical_partition():
    """Claim under attack: none — this is a regression check.

    Mutation: none (canonical DAG state).

    Structural expectation: I2 identity holds on the canonical
    residual partition {3, 16, 42} summing to 61.

    Observed behaviour required: I2 check's consistent flag True.
    """
    dag_reset()
    dag_put('C_total', 61, source='L_count', derivation='canonical')
    dag_put('d_eff', 102, source='L_self_exclusion', derivation='canonical')
    dag_put('n_baryon', 3, source='L_equip', derivation='canonical')
    dag_put('n_dark', 16, source='L_equip', derivation='canonical')
    dag_put('n_vacuum', 42, source='L_equip', derivation='canonical')

    acc = acc_SM()
    result = _identity_I2_gauge_cosmological(acc)

    assert result['consistent'], (
        f"I2 consistency failed on canonical partition: {result}")
    assert result['pi_F (K)'] == 61
    assert result['residuals_sum'] == 61
    assert result['Omega_m_plus_Omega_Lambda_closure']

    dag_reset()
    return {'passed': True, 'claim': 'I2 holds on canonical partition'}


# ──────────────────────────────────────────────────────────────────
# Test 3: I2 arithmetic fails loudly if residuals don't sum to K
# (so it's NOT a vacuous tautology)
# ──────────────────────────────────────────────────────────────────

def test_I2_fails_on_mutated_residuals_that_break_closure():
    """Claim under attack: I2 consistency is trivially true.

    Mutation: install residuals {3, 16, 43} — sum 62 != K=61. The
    closure arithmetic must fail.

    Structural reason: I2's residuals_sum_to_K check is an identity,
    not a tautology. Mutating the sum breaks the identity.

    Observed behaviour required: consistent flag False,
    residuals_sum_to_K False, Omega_m_plus_Omega_Lambda_closure False.

    This test proves I2 has non-vacuous content: it distinguishes
    valid partitions from invalid ones. A previous test ensured I2
    passes on the canonical partition; this test ensures I2 fails
    on a near-neighbour mutation, which together prove I2 is a real
    check, not a rubber stamp.
    """
    dag_reset()
    dag_put('C_total', 61, source='L_count', derivation='canonical')
    dag_put('d_eff', 102, source='L_self_exclusion', derivation='canonical')
    # Mutated: sums to 62, not 61
    dag_put('n_baryon', 3, source='test_smuggling', derivation='canonical')
    dag_put('n_dark', 16, source='test_smuggling', derivation='canonical')
    dag_put('n_vacuum', 43, source='test_smuggling',
            derivation='mutation: 42 -> 43, breaks residual closure')

    # ACC.__init__ itself has a partition-sum check that raises
    # ValueError on sum violations. That is stronger than I2 catching
    # the mismatch: smuggling is blocked at construction time, so the
    # invariant cannot be violated to even reach I2. Good.
    try:
        acc = acc_SM()
        # If ACC construction DID succeed, I2 should at least flag the
        # violation downstream.
        result = _identity_I2_gauge_cosmological(acc)
        assert not result['consistent'], (
            f"I2 consistency FALSELY passed on sum-violating partition — "
            f"means I2 is a rubber stamp, not a real check. {result}")
    except ValueError as e:
        # Expected: ACC refuses to construct on sum violations
        assert 'residuals sum' in str(e) and '!= K' in str(e), (
            f"ACC refused to build but for wrong reason: {e}")

    dag_reset()
    return {'passed': True,
            'claim': 'I2 has non-vacuous content; fails on sum-violating mutation'}


# ──────────────────────────────────────────────────────────────────
# Test 4: CANON-constant agreement with L_equip derivation
# ──────────────────────────────────────────────────────────────────

def test_CANON_residuals_agree_with_L_equip_derivation():
    """Claim under attack: the CANON constants (3, 16, 42) encode
    the answer the paper is supposed to derive.

    Mutation: run L_equip (which derives the residual partition from
    the matter-sector decomposition n_baryon=N_gen=3,
    n_dark=5*N_gen+1=16, n_vacuum=C_total-n_baryon-n_dark=42). Cross-
    check derived values against CANON constants.

    Structural expectation: the L_equip derivation produces the same
    three integers {3, 16, 42} that CANON hardcodes. If this test
    passes, the CANON constants are defensible as cached derived
    values, not as smuggled-in constants.

    Observed behaviour required: derived sectors match CANON exactly.

    Note: this test does NOT prove L_equip's derivation is itself
    free of upstream smuggling. Its 5*N_gen+1 identification of the
    dark sector depends on structural identifications that Phase 22
    will audit separately. This test only certifies that acc_SM()'s
    CANON constants match L_equip's output.
    """
    from apf.cosmology import check_L_equip

    dag_reset()
    # Bootstrap dependencies for L_equip
    dag_put('C_total', 61, source='L_count', derivation='canonical')
    dag_put('N_gen', 3, source='L_count', derivation='canonical')
    dag_put('d_eff', 102, source='L_self_exclusion', derivation='canonical')

    # Run L_equip — this populates n_baryon/n_dark/n_vacuum
    check_L_equip()

    assert dag_has('n_baryon'), "L_equip did not populate n_baryon"
    assert dag_has('n_dark'), "L_equip did not populate n_dark"
    assert dag_has('n_vacuum'), "L_equip did not populate n_vacuum"

    n_b = dag_get('n_baryon', consumer='test_smuggling')
    n_c = dag_get('n_dark', consumer='test_smuggling')
    n_v = dag_get('n_vacuum', consumer='test_smuggling')

    assert n_b == _CANON_OMEGA_B_NUM, (
        f"L_equip derives n_baryon={n_b} but CANON_OMEGA_B_NUM="
        f"{_CANON_OMEGA_B_NUM}. CANON constants have drifted from derivation.")
    assert n_c == _CANON_OMEGA_C_NUM, (
        f"L_equip derives n_dark={n_c} but CANON_OMEGA_C_NUM="
        f"{_CANON_OMEGA_C_NUM}. CANON constants have drifted.")
    assert n_v == _CANON_OMEGA_LAMBDA_NUM, (
        f"L_equip derives n_vacuum={n_v} but CANON_OMEGA_LAMBDA_NUM="
        f"{_CANON_OMEGA_LAMBDA_NUM}. CANON constants have drifted.")

    dag_reset()
    return {'passed': True,
            'claim': 'CANON residuals agree with L_equip derivation at v6.9'}


# ──────────────────────────────────────────────────────────────────
# Test runner
# ──────────────────────────────────────────────────────────────────

_TESTS = [
    test_acc_SM_reads_residuals_from_dag_when_populated,
    test_I2_passes_on_canonical_partition,
    test_I2_fails_on_mutated_residuals_that_break_closure,
    test_CANON_residuals_agree_with_L_equip_derivation,
]


def run_all_smuggling_tests():
    """Run every anti-smuggling test; return dict of results."""
    results = {}
    for test in _TESTS:
        name = test.__name__
        try:
            res = test()
            results[name] = {'status': 'PASS', **res}
        except AssertionError as e:
            results[name] = {'status': 'FAIL', 'error': str(e)}
        except Exception as e:
            results[name] = {'status': 'ERROR',
                             'error': f'{type(e).__name__}: {e}'}
    return results


if __name__ == '__main__':
    import sys
    results = run_all_smuggling_tests()
    print('=' * 72)
    print('Anti-smuggling test suite — Phase 22 pilot')
    print('=' * 72)
    n_pass = 0
    n_fail = 0
    for name, res in results.items():
        flag = res['status']
        print(f'[{flag:4s}] {name}')
        if 'claim' in res:
            print(f'         ↪ {res["claim"]}')
        if 'error' in res:
            print(f'         ↪ {res["error"]}')
        if flag == 'PASS':
            n_pass += 1
        else:
            n_fail += 1
    print('-' * 72)
    print(f'{n_pass} passed, {n_fail} failed')
    sys.exit(0 if n_fail == 0 else 1)
