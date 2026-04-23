---
type: concept
domain: apf
layer: 2-structure
created: 2026-04-14
updated: 2026-04-14
sources: []
---

# Gauge Uniqueness
## "The Standard Model Gauge Group is Forced"

**Core claim:** Given the [[Non-Closure Theorem|non-closure]] constraint and the requirement that the gauge structure must support 45 fermions (the necessary field content), the gauge group $SU(3) \times SU(2) \times U(1)_Y$ is the *unique* possibility. Not assumed; derived.

**Proof structure (L_gauge_template_uniqueness):**

1. **Non-closure rules out global symmetries** — L_NZ forbids any continuous automorphism that closes the Hilbert space. Classical gauge groups (e.g., $SU(5)$, $SO(10)$) would close the space under duality.

2. **Local gauge structure as cost-minimization** — The only admissible way to maintain symmetries is via local (spacetime-dependent) gauge fields. This is the only way to avoid closure.

3. **Representation constraints** — The 45 fermions must fit into irreps of the gauge group. The only factorization that works is:
   - Color: $SU(3)_c$ (octet representation for gluons; color triplets for quarks)
   - Isospin: $SU(2)_L$ (left-handed doublets; breaks electroweak)
   - Hypercharge: $U(1)_Y$ (quantizes fermion assignments)

4. **Why $N_c = 3$** (check_T_gauge) — The color number is forced by anomaly cancellation and the field content count.

5. **Why no higher-dimensional groups** — SU(4), SO(10), or E_6 would require additional fermions or imply global symmetries (violates L_NZ).

**Comparison:**
- Not a unification assumption (GUT theories *assume* larger groups and *derive* the standard model as a subgroup)
- APF *derives* the standard model group as the *only* admissible choice from first principles
- This is stronger: it explains why GUT unification is problematic from an admissibility perspective

**Status:** Proven in `gauge.py` (29 checks); formalized in [[Paper 2 - Structure]]

## See also
- [[Non-Closure Theorem]]
- [[Paper 2 - Structure]]
- [[Axiom A1]]
