---
type: synthesis
domain: apf
layer: cross-cutting
created: 2026-04-14
updated: 2026-04-17
sources: []
---

# Derivation Chain
## "The Logical Architecture of APF"

**Overview:** APF is built in seven nested layers, each building on the previous. This map shows the dependencies and key theorems at each step. Since the PLEC consolidation (April 2026), the Layer 0 foundation is presented as the *Principle of Least Enforcement Cost* — reality as the minimum-cost expression of distinction compatible with finite capacity — decomposed into four structurally necessary components: **A1** (finite capacity; upper bound), **MD** (positive cost floor), **A2** (argmin selection), **BW** (budget-window / non-degeneracy). Each is essential: Paper 1 v4.0-PLEC and its Technical Supplement §1 give countermodels showing no three of the four imply the fourth.

**Layer 0: Admissibility Geometry (PLEC four components)**
- [[Axiom A1]] — finite enforcement capacity (upper bound)
- MD — minimum cost floor (positive floor; makes "minimum" well-posed)
- A2 — argmin selection (picks realised structure from admissible set)
- BW — budget-window non-degeneracy (distinct-cost alternatives exist)
- L_epsilon_star — operability bound
- L_NZ, L_loc, L_cost, L_irr — lemmas on cost structure

**Layer I: Quantum Mechanics** (derives from L0)
- T2 → complex Hilbert space forced
- T3 → noncommutativity necessary
- T_Born → [[Born Rule]] derived
- T_Tsirelson → Bell inequality bound
- T_CPTP → quantum evolution as cost-preserving maps

**Layer II: Standard Model** (derives from I)
- [[Non-Closure Theorem]] (L_NZ) blocks classical gauge closure
- Theorem_R (R1/R2/R3) → representation constraints
- L_gauge_template_uniqueness → [[Gauge Uniqueness]]
- T_gauge → $SU(3) \times SU(2) \times U(1)$ forced, $N_c=3$
- T_field → 45 fermions necessary, no exotics
- L_count → exactly 61 effective degrees of freedom

**Layer III: Mass and Mixing** (derives from II)
- [[Gram Matrix]] — capacity structure determines spectrum
- L_multiplicative_amplitude → Feynman amplitude form
- T27c (x = 1/2) → binds exactly 3 generations
- L_mass_from_capacity → fermionic masses from cost topology
- T_CKM, T_PMNS → quark and lepton mixing angles

**Layer IV: Cosmology** (derives from I, II, III)
- L_equip → energy equipartition across sectors
- L_dark_budget → cosmological constant explained
- L_N_eff_prediction → neutrino family count
- L_equation_of_state → dark energy equation of state

**Layer V: Gravity and Spacetime** (derives from all)
- Internalization → Yang-Mills structure → Einstein equations (conjectured)
- Spacetime metric from cost topology
- Open: full derivation of Riemannian geometry from A1

**Interdependencies:**
- Layers I-II are logically tight; predictions follow deterministically
- Layers III-IV refine mass/cosmology predictions
- Layer V extends to gravity (partially complete)

**Implementation:** All 18 codebase modules, 294+ theorems, 349 checks

## See also
- [[Paper 13 - Minimal Admissibility Core]]
- [[What Physics Permits]]
- [[Predictions Catalog]]
