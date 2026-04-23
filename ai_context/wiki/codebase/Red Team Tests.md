---
type: entity
domain: apf
layer: cross-cutting
created: 2026-04-14
updated: 2026-04-14
sources: []
---

# Red Team Tests
## "Adversarial Checks Against APF"

**Purpose:** The red team module (`red_team.py`, 19 checks) systematically tests APF for contradictions, false positives, and hidden assumptions. If APF is wrong, these tests should catch it.

**Test categories:**

**1. Logical consistency (6 checks)**
- No contradiction between layers I-VI
- Gauge uniqueness is actually unique (no alternative solutions)
- Fermion count is determinate (not over-constrained)
- Mass ratios have unique solutions (Gram matrix is non-singular)
- No violations of admissibility axioms in derived structures

**2. Boundary and edge cases (5 checks)**
- Limiting case: $N \to 0$ (trivial system) behaves correctly
- Limiting case: $N \to \infty$ (infinite dimensions) is correctly forbidden
- Single-particle state ($d=1$ Hilbert space) is admissible baseline
- Two-particle entanglement satisfies Tsirelson bound
- Gauge group with $N_c \neq 3$ is provably inadmissible

**3. False positive detection (5 checks)**
- Naive cost bound gives wrong predictions (must refine to L_epsilon_star)
- Closure in subspace would be consistent (but full space is not closed)
- Fermion count is sensitive to representation choice (all give 45)
- No "accidental" agreement with experiment (predictions are structural)
- Dark energy prediction is unique (not degenerate with matter)

**4. Hidden assumption hunts (3 checks)**
- No appeal to symmetry principles (gauge structure is forced, not assumed)
- No circular reasoning in deriving Born rule (independent of quantum postulates)
- No anthropic selection in cosmological predictions
- Field content is not hand-fitted to data

**Results (v6.7):**
- All 19 checks **PASS**
- No contradictions found
- No hidden free parameters discovered
- No false positives identified
- Gauge uniqueness verified
- Fermion count over-determined (3 independent derivations all give 45)

**Failure modes (what would invalidate APF):**

If any of these happened, APF would be false:

1. Two different gauge groups both satisfy uniqueness conditions
2. The dark energy density had a different numerical value (would require new axiom)
3. Neutrino family count N_eff = 4.5 (was ruled out; still could appear in future data)
4. Fermion mass ratios showed large deviations (current: max 2.6%)
5. A viable classical unified group with closure (would violate L_NZ)

**How to extend:**

New red team tests can be added in `red_team.py` by:
1. Identifying a potential weakness
2. Formulating as a boolean check
3. Running against all 349 existing checks to ensure consistency
4. If test fails (returns False), a contradiction has been found

**Historical note:**

The red team module was developed *after* the core framework. It succeeded in finding several subtle errors in early versions (v3-v5) related to representation theory and cost accounting. Those errors have since been fixed. The fact that v6.7 passes all red team tests is a strong signal of internal consistency.

## See also
- [[APF Codebase]]
- [[Derivation Chain]]
- [[Open Problems]]
