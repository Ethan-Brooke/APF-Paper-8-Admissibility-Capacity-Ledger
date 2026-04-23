---
type: source
domain: related-physics
layer: cross-cutting
created: 2026-04-14
updated: 2026-04-14
sources: []
---

# Reconstruction Programs
## "Research Directions for Completing APF"

**Overview:** A collection of focused research programs aimed at resolving [[Open Problems]] and extending [[Admissible Physics Framework|APF]] into uncharted territory. Each program is a self-contained research direction with clear milestones and success criteria.

**RP1: Riemannian Geometry Derivation** (Priority S.1)

**Objective:** Derive Einstein field equations from [[Axiom A1]] via the admissibility cost topology.

**Current state:**
- Internalization gives Yang-Mills gauge bosons (formal, complete)
- Graviton as a composite object candidate (under development)
- Metric from cost-energy correspondence (conceptual)

**Missing pieces:**
1. Formal definition of cost in gravitational sector
2. Proof that Einstein equations minimize cost subject to admissibility
3. Connection to Newton's constant from admissibility parameters

**Success criterion:** Derive $R_{\mu\nu} - \tfrac{1}{2}R g_{\mu\nu} + \Lambda g_{\mu\nu} = \kappa^2 T_{\mu\nu}$ from A1 alone.

**Timeline:** 12-18 months (requires new check functions in gravity.py)

---

**RP2: Dark Energy Equation of State** (Priority A.2)

**Objective:** Extend dark energy sector to accommodate both $w = -1$ and dynamical $w(z)$ scenarios.

**Current tension:**
- APF predicts $w = -1.00$ (cosmological constant)
- DESI DR2 suggests $w \approx -0.9$ (possible evolving dark energy)

**Two pathways:**

*Path A:* DESI systematic error (most likely)
- Final DESI data contradicts current hints
- APF prediction stands
- Timeline: Q2 2025 when final data released

*Path B:* Dynamical dark energy (alternative)
- Coupled scalar field in admissibility framework
- Potential $V(\phi)$ derived from cost principles
- Attractor solution reaches $w \approx -0.95$ at late times
- Timeline: 6 months if pursued

**Success criterion:** Either DESI confirms $w = -1.00$, or dynamical model matches data to 0.5%

---

**RP3: Matter-Antimatter Asymmetry Mechanism** (Priority A.4)

**Objective:** Formalize CP violation in the early-universe admissibility sector.

**Current state:**
- η_B prediction to 0.49% accuracy
- Mechanism sketched in Paper 7; not fully implemented

**Missing implementation:**
1. Detailed CP-violating operators in early admissibility sector
2. Numerical integration of asymmetry evolution
3. Connection to CKM and PMNS phases

**Candidate mechanism:**
- Kaon oscillation analog in admissibility Hilbert space
- Cost asymmetry in forward vs. reverse baryon number violation
- Predicts specific baryon decay modes

**Success criterion:** Derive η_B with <1% error; predict observable rare decay channels

**Timeline:** 6-9 months

---

**RP4: Dark Matter Particle Model** (Priority A.5)

**Objective:** Identify the dark matter particle from admissibility constraints.

**Candidates under consideration:**

1. **Sterile neutrino** (favored by APF)
   - Mass range: 0.1–1 keV (from N_eff and cosmological constraints)
   - Coupling: gravity-only (avoids X-ray constraints)
   - Production: non-thermal in early universe
   - Predicted signal: None (dark)

2. **Axion** (sub-dominant)
   - Can be generated from admissibility; less constrained than sterile neutrino
   - ADMX/CAST experiments could detect

3. **Composite state** (exploratory)
   - Bound state of admissible fermions with non-standard interactions
   - Novel, testable signatures

**Success criterion:** Identify dominant dark matter species; make testable prediction (detection rate, annihilation cross-section, or decay mode)

**Timeline:** 12 months

---

**RP5: Quantum Gravity Unification** (Advanced)

**Objective:** Merge quantum field theory and general relativity via admissibility.

**Dependencies:**
- Completes RP1 (Riemannian geometry)
- Requires resolution of RP2 (dark energy)

**Approach:**
- Admissibility constraints at Planck scale
- Cost budget between matter, gauge, and gravitational sectors
- No separate gravity "quantization" needed; gravity emerges from cost principle

**Success criterion:** Consistent quantum field theory on curved spacetime; predictions for Planck-scale physics

**Timeline:** 24–36 months (post-RP1)

---

**RP6: Inflation Mechanism** (Priority B.7)

**Objective:** Derive inflation from early-universe admissibility cost profile.

**Current state:**
- Conceptual sketch in Paper 7
- No detailed slow-roll predictions

**Mechanism:**
- Early universe has higher enforcement cost per degree of freedom
- This creates a "cost barrier" analog to inflaton potential
- Inflation ends when cost structure transitions (epoch change)

**Success criterion:** Predict tensor-to-scalar ratio r, spectral index n_s, running α_s within 5% of observations

**Timeline:** 9-12 months

---

**Schedule (v6.7 → v8.0):**

```
Now (2026-Q2)       RP2 monitoring (DESI), RP4 start
2026-Q3             RP3 formalization, RP1 preliminary
2026-Q4             RP1 formal proofs, RP4 completion
2027-Q1             RP2 resolution (DESI final), RP5 planning
2027-Q2–Q4          RP1 completion, RP5 initial development
2028+               RP5 quantum gravity, RP6 if time permits
```

## See also
- [[Open Problems]]
- [[Paper 7 - Action]]
- [[APF Codebase]]
