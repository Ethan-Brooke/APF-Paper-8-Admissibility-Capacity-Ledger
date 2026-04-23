---
type: concept
domain: apf
layer: 1-spine
created: 2026-04-14
updated: 2026-04-17
sources: []
---

# Axiom A1
## "Any physical distinction that can be maintained must be maintainable by a process with finite operational cost."

> **PLEC note (2026-04-17):** As of the PLEC consolidation (Papers 1 v4.0-PLEC and 2 v5.3-PLEC), A1 is one of **four structurally necessary components** of the *Principle of Least Enforcement Cost* (PLEC): A1 (finite capacity; upper bound), **MD** (positive cost floor), **A2** (argmin selection), **BW** (budget-window / non-degeneracy). PLEC is the canonical variational formulation of the enforceability-of-distinction foundation; it does not replace the locked claim that physical meaning is grounded in enforceable distinction. Each component has an explicit countermodel / essentiality proof in the Paper 1 Technical Supplement showing it is not implied by the other three.

> **Two-axiom history (2026-04-16, superseded by PLEC framing):** Before the PLEC consolidation, A1 and A2 were presented as two named axioms. Post-PLEC, A1 and A2 are both components of a single variational principle. The "A1 admits / A2 selects" distinction remains correct as a description of which component does which work: A1 defines the admissible set $\mathcal{A}(\rho, R)$ of configurations satisfying $\Sigma\varepsilon \le C$; A2 picks the min-cost member as the realised structure. MD is the floor that makes "minimum" well-posed; BW is the non-degeneracy ensuring distinct-cost alternatives exist.

**Statement:** Let $D$ be a physical distinction (a difference between two possible states of a system). If $D$ can be maintained indefinitely by some operational process, then that process must have finite cost (in the sense of resources, energy, or work required per unit time).

**Interpretation:**
- "Physical distinction" — a state difference that (a) can be prepared, (b) can be measured, (c) has physical consequences
- "Maintained" — the distinction persists or can be repeatedly enforced
- "Operational process" — an implementable, repeatable algorithm or natural law
- "Finite cost" — the resource requirement (energy, work, logical operations) is bounded per cycle

**Consequences:**

From PLEC (A1 + MD + A2 + BW), APF derives:
1. Complex Hilbert space structure (not real or quaternionic) — A1
2. The [[Born Rule]] probability formula — A1
3. $SU(3) \times SU(2) \times U(1)$ gauge group uniquely — A1 admits the gauge templates; A2 selects $N_c = 3$ as the min-cost member
4. 45 fermions in 3 families — A1 + A2 (1-of-4680 enumeration; A2 picks the min-cost survivor)
5. Exact mass ratios and mixing angles
6. Dark energy and cosmological density fractions
7. 47+ quantitative predictions with zero free parameters
8. $\theta_{QCD} = 0$ — A2 selects $\theta = 0$ over $\theta \ne 0$ (both admissible under A1; A2 picks zero-maintenance-cost value)

**Comparison to other axioms:**
- Not a variational principle (doesn't assume least-action)
- Not tied to symmetry (does not assume gauge invariance)
- Precedes and forces gauge structure (gates the latter)
- Ground-level: applies to any system maintaining distinguishability

**Key lemmas (from [[Paper 1 - Spine]]):**
- L_epsilon_star — operability bound on cost
- L_NZ — non-closure of state spaces
- L_loc — locality from finiteness
- L_cost — cost scaling with ensemble size
- L_irr — irreversibility from cost asymmetry

**Status:** Foundational; fully formalized in `core.py` (48 checks)

## See also
- [[Enforcement Budget]]
- [[Paper 1 - Spine]]
- [[Non-Closure Theorem]]
