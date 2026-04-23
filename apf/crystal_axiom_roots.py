"""APF v6.9 — Phase 13.1 anchor: PLEC component → source-most bank check.

This module is a pure-data constant. It encodes the canonical mapping from
each of the four PLEC ontological/structural commitments referenced by the
Enforcement Crystal (Paper 14, v2 corrected) and the Admissibility-Capacity
Ledger (apf.unification) to its *source-most* bank-registered check in the
v6.9 codebase.

"Source-most" means: the bank entry that *first* introduces the commitment
into the framework (as an axiom, lemma, or executable witness) — the entry
that all downstream uses ultimately reduce to. Downstream papers and modules
(Papers 4–7, apf/unification.py, apf/unification_three_levels.py,
apf/crystal.py) reference these checks transitively via dependency strings.

PLEC components
---------------
A1  Finite Enforcement Capacity (Axiom):
    sum_{d in D(rho,R)} epsilon(d) <= C(R) < infinity.
    Source-most: check_A1 (apf/core.py).
    Status: postulated (the axiom itself).

A2  Minimum Cost Selection (Axiom, formerly implicit):
    G_realized = argmin_{q in A_Gamma} K[q].
    Source-most: check_Regime_R (apf/plec.py).
    Status: [P] — well-posedness of the selector under R1..R4
    (smooth, locally additive, connected, unsaturated).
    Logically independent of A1 + MD (countermodels in Paper 1
    Supplement v2 §1).

MD  Uniform Per-Test Cost Floor (Structural Assumption):
    epsilon(d) >= mu^* * n(d) for some mu^* > 0,
    where n(d) is the per-test multiplicity carried by d.
    Source-most: check_L_epsilon_star (apf/core.py).
    Status: [P] — the n(d) = 1 case is proved by compactness;
    the multi-test extension is a uniform-bound restatement.
    Without MD, L_eps*, T_form, L_cost, and the Wedderburn-Artin
    chain all fail (Paper 1 Supplement v2 §1).

BW  Budget-Window Richness (Richness Premise):
    The cost spectrum at a single interface contains a witness
    triple (d1, d2, d3) such that {d1, d3} is admissible but
    {d2, d3} is not — i.e. the budget window discriminates
    second-distinction admissibility.
    Source-most: check_worked_example (apf/core.py).
    Status: [P] — explicit witness with C = 5,
    eps_1 = 2, eps_2 = 3, eps_3 = 5/2 demonstrating the
    triple condition.
    Holds generically whenever the cost spectrum carries
    >= 3 distinct values (Paper 1 Supplement v2 §1).

Use
---
This module exposes a single constant ``PLEC_AXIOM_ROOTS`` plus a thin
lookup helper ``axiom_root(name)``. Downstream modules (in particular
``apf/crystal.py``, planned for Phase 13 Stage II) consume these to
build the Enforcement Crystal's axiom-root vertex set.

This module is *not* itself bank-registered — it carries no executable
check. The four checks it points at are already individually registered
(check_A1 in core, check_Regime_R in plec, check_L_epsilon_star in core,
check_worked_example in core), so verify_all already exercises them.

References
----------
- Paper 14 (Enforcement Crystal v2 corrected, Zenodo 18615555)
- Paper 1 Technical Supplement v2 §1 (PLEC essentiality proofs;
  countermodels for A1 vs A2 vs MD logical independence)
- apf/unification.py header (cross_refs list at line 816 lists
  the same four anchors A1, MD, BW, A2 as ontological/structural
  commitments)
- Phase 13 Stage I task 13.1 in
  ``__APF Library/APF Reference Docs/Reference - APF Paper Update Work Plan v2.md``
"""

from collections import OrderedDict


# =============================================================================
# Canonical mapping: PLEC component → source-most bank check
# =============================================================================

PLEC_AXIOM_ROOTS = OrderedDict([
    ("A1", {
        "name":        "Finite Enforcement Capacity",
        "kind":        "axiom",
        "statement":   "sum_{d in D(rho,R)} epsilon(d) <= C(R) < infinity",
        "check":       "A1",
        "check_func":  "check_A1",
        "module":      "apf.core",
        "module_path": "apf/core.py",
        "epistemic":   "axiom",
        "essentiality_ref": "Paper 1 Supplement v2 §1 box 'A1'",
    }),
    ("A2", {
        "name":        "Minimum Cost Selection",
        "kind":        "axiom",
        "statement":   "G_realized = argmin_{q in A_Gamma} K[q] (PLEC selector)",
        "check":       "Regime_R",
        "check_func":  "check_Regime_R",
        "module":      "apf.plec",
        "module_path": "apf/plec.py",
        "epistemic":   "P",  # selector well-posedness is proved on R1..R4
        "essentiality_ref": "Paper 1 Supplement v2 §1 box 'A2'",
    }),
    ("MD", {
        "name":        "Uniform Per-Test Cost Floor",
        "kind":        "structural_assumption",
        "statement":   "epsilon(d) >= mu^* * n(d), mu^* > 0",
        "check":       "L_epsilon*",
        "check_func":  "check_L_epsilon_star",
        "module":      "apf.core",
        "module_path": "apf/core.py",
        "epistemic":   "P",  # n(d)=1 case proved by compactness
        "essentiality_ref": "Paper 1 Supplement v2 §1 box 'MD'",
    }),
    ("BW", {
        "name":        "Budget-Window Richness",
        "kind":        "richness_premise",
        "statement":   "exists triple (d1,d2,d3) s.t. {d1,d3} admissible, {d2,d3} not",
        "check":       "worked_example",
        "check_func":  "check_worked_example",
        "module":      "apf.core",
        "module_path": "apf/core.py",
        "epistemic":   "P",  # explicit witness in worked example
        "essentiality_ref": "Paper 1 Supplement v2 §1 box 'BW'",
    }),
])


# =============================================================================
# Lookup helper
# =============================================================================

def axiom_root(name):
    """Return the PLEC anchor record for ``name`` (one of 'A1', 'A2', 'MD', 'BW').

    Raises KeyError on an unknown name.
    """
    if name not in PLEC_AXIOM_ROOTS:
        raise KeyError(
            f"Unknown PLEC component {name!r}. "
            f"Expected one of: {list(PLEC_AXIOM_ROOTS)}"
        )
    return PLEC_AXIOM_ROOTS[name]


def all_roots():
    """Return the four-element list of (name, record) pairs in canonical order."""
    return list(PLEC_AXIOM_ROOTS.items())
