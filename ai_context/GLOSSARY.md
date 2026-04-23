# GLOSSARY.md — APF Vocabulary

Terms are grouped by category. Each entry is short; follow the cross-references in `ai_context/FRAMEWORK_OVERVIEW.md` and the main paper for depth.

---

## Axioms and structural components

**A1 (Finite Information Capacity).** The project's single axiom. Every interface of physical distinction has finite information capacity. Formally: $\sum_{d \in G} \varepsilon(d) \leq C(\rho, \Gamma)$ for some finite $C$. Upper bound on the admissible set.

**MD (Minimum Distinction cost).** Structural companion to A1. Every distinction has strictly positive enforcement cost: $\varepsilon(d) \geq \mu^* > 0$. Lower bound on the cost spectrum.

**A2 (Argmin selection).** Realized physical structure is the least-enforcement-cost member of the admissible class. The "realization" component of PLEC.

**BW (Non-degeneracy).** The cost spectrum is non-degenerate — no two distinct distinctions have identical cost functionals. Spectrum discipline.

**A4 (Irreversibility).** Derived consequence. Cross-interface correlations commit capacity that no local observer can recover. Equivalent to: the arrow of time emerges from locality + finite capacity. See `L_irr` below.

**A9.1–A9.5.** The five geometric prerequisites for Lovelock uniqueness in $d = 4$ (locality, covariance, conservation, second-order, propagation). Derived in v6.9 via `check_A9_closure` in `apf/gravity.py`.

## The PLEC (Principle of Least Enforcement Cost)

**PLEC.** The four-component variational principle: A1 (upper bound) + MD (lower bound) + A2 (argmin selection) + BW (non-degeneracy). Not a new axiom; a canonical statement of the variational structure of the framework.

**Admissible class $\mathcal{A}_\Gamma$.** The set of configurations physically licensable in context $\Gamma$ under A1 + MD + BW.

**Enforcement-cost functional $E_\Gamma$.** Maps admissible configurations to their total capacity burden: $E_\Gamma: \mathcal{A}_\Gamma \to [0, \infty)$.

**PLEC selection $X^*$.** The realized configuration: $X^* \in \arg\min_{X \in \mathcal{A}_\Gamma} E_\Gamma(X)$.

**Regime R.** The regime in which PLEC selection is well-posed: R1 smooth cost, R2 local additivity, R3 connected path space, R4 no saturation. Formalized as `check_Regime_R` in `apf/plec.py` (v6.9).

## Regime-exit taxonomy (five structurally distinct exits)

Formalized as `check_Regime_exit_Type_I` through `check_Regime_exit_Type_V` in `apf/plec.py` (v6.9).

**Type I — collapse of admissible variation.** Saturation. Admissible neighborhood shrinks to zero measure; PLEC selection becomes trivially unique because no variation remains.

**Type II — minimizer nonuniqueness.** Branching. Admissible class is nonempty but the least-cost selector has multiple inequivalent minimizers (up to the relevant equivalence relation).

**Type III — change of admissible class.** Record locking. The relevant admissible class itself changes (e.g., coherent $M_{\mathrm{sys}}$ → record-locked $M_{\mathrm{sys}} \otimes Z_R$). Transition is not a breakdown internal to one class but a transfer to a different domain.

**Type IV — loss of smooth or local structure.** Singularities, Planck-scale discreteness, topology change. Admissibility is intact; regularity assumptions (smoothness, tangent-space structure, chartability) fail.

**Type V — pure representational redundancy.** Gauge freedom, coordinate ambiguity. Admissibility and the realized minimizer are intact; only the descriptive coding is non-unique. NOT physical branching (distinct from Type II).

## Capacity and partition

**$C_{\mathrm{total}} = 61$.** The total type count of the APF matching. Fixed integer. Not tunable.

**The 61-partition.** $61 = 42 + 3 + 16 = 42 + 19$. Vacuum types (42) + baryonic matter types (3) + cold-dark-matter types (16). Matter sub-total is 19. The partition is topological (type-counting), not dynamical. `check_L_saturation_partition`.

**$C_{\mathrm{ext}}$ / $C_{\mathrm{int}}$.** External (geometric) capacity vs. internal (gauge-addressable) capacity. $C_{\mathrm{int}} = 12 = \dim(SU(3) \times SU(2) \times U(1))$. $C_{\mathrm{ext}}$ supports geometric degrees of freedom; in Paper 6, identified as the geometric-sector candidate for dark matter phenomenology (tagged [P_structural] + [Pred]).

**Capacity residual (Λ).** T11: the cosmological constant is the global capacity residual after all local enforcement commitments, $\Lambda \sim (C_{\mathrm{total}} - C_{\mathrm{used}})/V \sim 1/R_H^2$.

## Epistemic tags

**[P] — Internally derived.** Result follows from APF axioms plus explicitly stated regime assumptions, without external classification or uniqueness theorems.

**[P] + [I] — Derived with import.** APF prepares an admissible class; an external mathematical theorem (Solèr, Gleason, HKM, Malament, Lovelock, etc.) closes the classification or uniqueness. The framework is responsible for preparing the class; the import closes it.

**[P_structural] — Structural identification.** The framework identifies an already-derived structure with a familiar physical role. Weaker than a full derivation; often definition-level once prerequisites are in place.

**[Pred] — Phenomenological prediction.** An empirical claim downstream of the prior layers, to be evaluated against data. Scaling laws, sector identifications, exclusion claims.

**[P + lattice]** — Proved using lattice QCD as an external input. Similar status to [P] + [I] but specifically calling out lattice-derived values.

**[C] — Conjecture.** Unproved. Rare in the current codebase (1 of 342 in v6.9). Appears only where no proof has been constructed yet.

## Core lemmas / theorems referenced throughout

**L_irr.** Irreversibility from cross-interface correlations. Cross-interface correlations commit capacity at both interfaces; no local observer can recover it. Monotonicity (L3) holds at each interface. **Note:** the v6.8/v6.9 framing is "locally unrecoverable" — NOT the older v4.x framing "reversal requires infinite cost" (which contradicts L3 and has been rejected). `check_L_irr` in `apf/core.py`.

**L_nc (non-closure).** Composition of interfaces can fail — not every pair of systems is jointly enforceable. The root of superadditive cost.

**L_col.** Minimality / colour-type constraint arising in the gauge derivation chain.

**L_saturation_partition.** The 61-partition is topological: type-classification by gauge-addressability ($Q_1$) and confinement ($Q_2$) is determined by discrete type labels, not by capacity amounts. Saturation-independent.

**L_anomaly_free.** The 61 types must coexist simultaneously. Anomaly cancellation (7 independent conditions on hypercharges) fails if any type is removed.

**L_equation_of_state.** $w = -1$ exactly at all cosmic epochs. Four escape routes enumerated and each blocked by a [P] theorem.

**L_equip.** At any saturation $s > s_{\mathrm{crit}}$, max-entropy distributes available capacity uniformly over the 61 types. Surplus-independent density fractions.

**L_Gleason_finite.** APF-internal verification of Gleason's theorem in finite dimension. Full Gleason remains an external [I] import.

**T_gauge.** Selection of $SU(3) \times SU(2) \times U(1)$ as the unique admissible gauge group.

**T_field.** The 45-chiral-fermion content (1 of 4680 anomaly-free sets survives).

**T4F.** $N_{\mathrm{gen}} = 3$ exactly, from capacity saturation.

**T7B.** Metric tensor from non-factorization + polarization identity. Quadratic non-factorizing response at shared interfaces IS a metric tensor.

**T8.** $d = 4$ as the minimal admissible spacetime dimension. $d \leq 3$ hard-excluded; $d \geq 5$ disfavored.

**T9 / L3-µ.** Record-locking channel formalization. For $k$ enforcement operations applied in all $k!$ orderings, the resulting CP maps are pairwise distinct. Yields $|H_k / \equiv| = k!$ inequivalent histories. Algebraic model for measurement.

**T9_grav.** Einstein equations from Lovelock closure on A9.1–A9.5 prerequisites. Unique in $d = 4$.

**T10.** $\kappa \sim 1/C_*$. Gravitational coupling as capacity scaling.

**T11.** Cosmological constant as global capacity residual. $\Lambda \sim (C_{\mathrm{total}} - C_{\mathrm{used}})/V$.

**T_horizon_reciprocity.** Horizon dissolves reciprocity: L_irr applied at each crossing means the bulk-symmetric pairing fails at horizons.

**T_Bek.** Bekenstein area-law bound. Entropy on an interface scales with area, not volume.

## Imported closures (the [I] imports)

APF prepares the admissible class; these external theorems close the uniqueness or classification.

**Solèr's theorem (1995).** Infinite-dimensional orthomodular spaces over a division algebra reduce to $\mathbb{R}, \mathbb{C}, \mathbb{H}$. Continuous one-parameter reversible families in 2D subspaces select $\mathbb{C}$.

**Gleason's theorem (1957).** In a complex Hilbert space of dimension $\geq 3$, the unique probability measure satisfying refinement additivity and admissibility invariance is $p(E|\rho) = \text{Tr}(\rho E)$. APF-internal preparation: `check_L_Gleason_finite`.

**Hawking–King–McCarthy theorem (1976).** A strict partial order on events satisfying density conditions determines conformal structure. APF preparation: A4 (irreversibility) gives the partial order.

**Malament's theorem (1977).** Conformal structure + volume normalization (from Radon–Nikodym uniqueness) fixes Lorentzian signature.

**Lovelock's theorem (1971).** In $d = 4$, the unique symmetric divergence-free tensor depending linearly on second derivatives of the metric is the Einstein tensor plus cosmological term. APF preparation: A9.1–A9.5 via `check_A9_closure`.

**Kolmogorov consistency (1933).** Continuum structure from finite-additivity + consistency. APF use: smooth-manifold emergence from discrete admissibility.

**Ostrogradsky instability.** Higher-than-second-derivative Lagrangians are unstable. APF use: justifies A9.4 (second-order only).

## Recurring symbols

| Symbol | Meaning |
|---|---|
| $\Gamma$ | Interface / resource context |
| $\rho$ | State |
| $\varepsilon(d)$ | Enforcement cost of distinction $d$ |
| $\mu^*$ | Minimum distinction cost (from MD) |
| $C$ or $C(\rho, \Gamma)$ | Interface capacity |
| $\mathcal{A}_\Gamma$ | Admissible class |
| $E_\Gamma$ | Enforcement-cost functional |
| $\Omega_X$ | Sector density fraction (X ∈ {Λ, m, b, cdm}) |
| $\epsilon_*$ or $\varepsilon_*$ | Quantum of cost at the de Sitter endpoint |
| $d_{\mathrm{eff}}$ | Effective internal dimension ($= 102$) |
| $\Lambda_{\mathrm{eff}}(k)$ | Effective cosmological constant at saturation level $k$ |
| $s$ | Saturation parameter |
| $s_{\mathrm{crit}}$ | Critical saturation threshold ($= 1/d_{\mathrm{eff}} = 1/102$) |
| $M_{\mathrm{sys}}$ | Coherent system algebra |
| $Z_R$ | Classical log algebra (record-locked regime) |

## Usage-warning terms

These terms have specific meanings on APF that differ from general usage. Getting them wrong misleads.

**"Quantum structure."** In APF: the specific Hilbert-space + Born-rule + CPTP representation that emerges as the realized least-cost bookkeeping in the coherent admissible regime. Not synonymous with "quantum mechanics" in the general sense; APF's claim is conditional on the coherent regime and exits at record-locking.

**"Measurement."** In APF: Type III regime exit (change of admissible class, coherent → record-locked). Not a collapse axiom; a formal class change via `check_T9`.

**"Dark energy."** In APF: the cosmological constant, with $w = -1$ exactly. Not a dynamical field. The four escape routes to $w \neq -1$ are all blocked by [P] theorems.

**"Dark matter."** In APF: the $C_{\mathrm{ext}}$ geometric capacity sector. Gravitationally active, gauge-inert. Phenomenological identification tagged [P_structural + Pred]; not a particle species.

**"Irreversibility."** In APF: locally unrecoverable cross-interface correlations (v6.8+ framing). NOT "infinite-cost reversal" (v4.x framing, rejected because it contradicts L3 monotonicity).

---

*Generated by the `create-repo` skill. To update, edit `templates/GLOSSARY.md.template` in the skill folder.*
