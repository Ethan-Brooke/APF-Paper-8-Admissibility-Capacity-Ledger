# OPEN_PROBLEMS.md — Known-Open APF Research Problems

This file lists the framework-level open problems as of the release date of this repo. The purpose is specifically to prevent AI agents (or external proposers) from naively re-solving known-open problems with shallow mechanisms.

If you are proposing new APF work, check here first. A proposal that claims to "resolve" a problem listed below without addressing the verdict on each escape route (where catalogued) is almost certainly reinventing the specific failure mode the skill `audit-response` exists to prevent.

---

## 1. The $H_0$ tension / dynamical dark energy

### Current APF prediction
$H_0 = 67.76$ km/s/Mpc (derived from $\Omega_m = 19/61$ and Planck's measured $\Omega_m h^2 = 0.143$, zero free parameters). $w = -1$ exactly at all cosmic epochs (`check_L_equation_of_state` [P]).

### Current observational situation
- **H0DN 2026** (Casertano et al., A&A 708 A166, published 10 April 2026): $H_0 = 73.50 \pm 0.81$ km/s/Mpc. Tension with APF prediction: $7.09\sigma$.
- **DESI DR2**: prefers $w_0 \approx -0.75$, $w_a \approx -0.99$. Tension with $w = -1$: $\sim 1.5$–$2\sigma$.

### Escape routes and verdicts
The $H_0$ tension can in principle be resolved by one of five mechanisms. Each has an APF verdict.

| Route | Mechanism | APF blocker | Empirical viability |
|---|---|---|---|
| I | Pre-recombination $r_d$ cut $\sim 7.8\%$ (early dark energy, extra $N_{\mathrm{eff}}$) | BLOCKED by L_equation_of_state, L_anomaly_free | EDE is mainstream candidate |
| II | Late-time dynamical DE ($\Omega_\Lambda(a)$ growth $\sim 12\%$) | BLOCKED by L_saturation_partition | DESI prefers wrong sign |
| III | Local distance-ladder systematics | none | H0DN 2026 substantially rules out |
| IV | Spatial curvature $\Omega_k \sim -0.05$ | none | Planck $\|\Omega_k\| < 0.02$, insufficient |
| V | Local underdensity $\sim 25\%$ at Hubble-flow scale | none | Observed $\sim 5$–$15\%$, short by $\sim 2\times$ |

No route is both APF-admissible and currently empirically sufficient. Route V is the only unblocked-and-not-yet-ruled-out APF-side path, and it's insufficient by a factor of $\sim 2$ against current observations.

### Falsifier statement
If H0DN's result stands as the true local expansion rate and is not Routes III/IV/V, then at least one of these [P] theorems must be demoted in a future codebase version:
- `check_L_equation_of_state` (would have to admit $w \neq -1$ in some regime)
- `check_L_saturation_partition` (would have to admit non-topological partition evolution)
- `check_L_anomaly_free` (would have to admit late-time type decoupling)

### What has been rejected
An external proposal (Grok, April 2026) attempted to resolve this via time-locked vacuum ($L(a) = L_{\max}(1 - a^{-\beta})$) producing phantom $w(z) \approx -1.5$. **Rejected** because:
1. Has two free parameters presented as zero ($\beta$, $L_{\max}$).
2. Moves in the wrong direction vs. DESI.
3. Implicitly demotes L_saturation_partition without acknowledging it.
4. Scaling $\rho_{\mathrm{DE}} \propto \Omega_\Lambda^2$ asserted, not derived.
5. $L(a)$ functional form ad hoc.
6. No end-to-end numerical recovery of $H_0 = 73.5$.

Full rejection memo in the canonical repo: `APF Reference Docs/Reference - H0 Tension Mathematical Analysis (2026-04-18).md`.

### Open work
- Route V (local inhomogeneity) is open; observational programs (KBC, Whitbourn-Shanks, Hoscheit-Barger) may deliver sharper numbers.
- Any proposed APF-side mechanism must (a) derive the functional form rather than assume it, (b) explicitly argue for which [P] theorem demotion is required if any, (c) check the sign/magnitude against DESI.

---

## 2. Absolute mass scales for $m_t$ and $m_b$

### What's derived
Structural mass ratios and hierarchy: the mass matrix texture, the seesaw structure, and the relative scales match the Standard Model observations within a few percent.

### What's open
The absolute scale for the top-quark mass and bottom-quark mass ($m_t, m_b$ in GeV). Currently: APF ratios fit but the overall scale ($v \sim 246$ GeV connection) requires a derivation of $\sigma$ (the dimensional scale parameter) that is not yet in hand.

### What would close it
A derivation of $\sigma$ from capacity / PLEC, analogous to how other scales (Planck scale, Higgs VEV) would emerge from the framework at appropriate regime boundaries.

### What would invalidate a proposed fix
Any proposed $m_t$ / $m_b$ scale derivation that leaves existing mass-ratio checks unchanged and produces a scale without touching $\sigma$ is probably sneaking a free parameter in under a different name.

---

## 3. $m_c$ at 2.6%

### What's derived
Ratios of charm-quark mass to other masses; structural consistency within the second generation.

### What's open
A persistent 2.6% offset between APF prediction and observation for $m_c$.

### Current status
Believed irreducible within the Schur structural limit: the 2.6% is the granularity of the current structural derivation. Improvement would require a deeper derivation of the Schur structure.

### Implication
Don't try to improve $m_c$ alone without also improving the Schur structural derivation.

---

## 4. Dark matter particle identity

### What's derived
- T12: derivation of the existence of a sector that gravitates but carries no gauge charge.
- Gravitational properties (rotation curves, lensing-compatible).
- `C_ext` identification as the geometric-sector candidate (tagged [P_structural] + [Pred]; see Paper 6 §13).

### What's open
Particle identity. APF does not identify dark matter as any specific Standard Model or BSM particle. The framework-internal candidate is "geometric capacity sector," which is a different kind of answer than "a particle species."

### What would close it
A bridge theorem: given the APF external sector $C_{\mathrm{ext}}$, derive the effective source term it contributes in the geometric response equations and prove that its admissible distributions reproduce at least one nontrivial dark-matter data class (rotation curves, lensing, structure formation) without introducing particle degrees of freedom.

### What would invalidate a proposed fix
Any proposed dark-matter identification that requires the 61-partition to change (e.g., a new type outside the 3 + 16 + 42 split) will break `L_anomaly_free` and is not admissible.

---

## 5. PLEC existence / uniqueness theorem

### What's stated
The four-component variational principle (A1 + MD + A2 + BW) is the canonical framing of APF's selection structure.

### What's open
A rigorous existence and uniqueness theorem for the PLEC selector: given $\mathcal{A}_\Gamma$, the enforcement-cost functional $E_\Gamma$, and the equivalence relation $\sim$, under what conditions does $\arg\min E_\Gamma$ exist and is unique up to $\sim$?

### What would close it
A theorem analogous to the direct method of the calculus of variations, with APF-specific coercivity / lower-semicontinuity analogues.

### Paper 6 status
Listed as Project C.1 in Paper 6 v2.0-PLEC Appendix C. `check_Regime_R` (v6.9) formalizes the regime in which PLEC is well-posed but does not fully prove existence + uniqueness in general.

---

## 6. Formal regime-exit theorem

### What's stated
The five-type regime-exit taxonomy (Types I–V) is formalized as bank checks in `apf/plec.py` (v6.9). Each type has an executable witness.

### What's open
A unified theorem tying each regime exit to a precise PLEC-hypothesis breakdown. Currently the taxonomy is structural: each Type is demonstrated on a specific witness. A general theorem of the form "if PLEC hypothesis $H_i$ fails, the corresponding exit is of Type $i$" would upgrade the taxonomy from vocabulary to machinery.

### Paper 6 status
Listed as Project C.4 in Paper 6 v2.0-PLEC Appendix C.

---

## 7. Dark-energy observational confrontation at DR3+

### What's stated
`check_L_DESI_response` quantifies the current $w_0$/$w_a$ tension with DESI DR2 at 1.5–2σ.

### What's open
Response to DESI DR3 (expected ~2027). If DR3 confirms $w_0 \neq -1$ at $\geq 5\sigma$, framework demotion is triggered.

### What will have been required
- Demotion of `check_L_equation_of_state` from [P] to regime-conditional.
- Rederivation with the added regime conditions.
- Revision of Paper 6 §11 ($\Lambda$ as capacity residual).

### Don't try to pre-empt this
Dynamical-DE proposals before DR3 results land are not admitted. The framework's current position is stated; revisions wait for the data.

---

## 8. Neutrino sector — structural ratios vs. phase derivation

### What's derived
- Neutrino hierarchy structure.
- Seesaw factorization.
- PMNS structural constraints matching observed magnitudes.

### What's open
- The CP phase $\delta_{\mathrm{PMNS}}$ derivation (session `delta_pmns` in progress).
- The absolute neutrino-mass scale.

### Paper status
In-session module `session_delta_pmns` tracks progress. Not yet in a stable paper.

---

## 9. Quantum-gravity / Planck-scale unification

### What's stated
APF derives GR in $d = 4$ (Paper 6) and QM in the coherent unsaturated regime (Paper 5). Both exit at their respective regime boundaries (horizons, Planck-scale discreteness, measurement, etc.).

### What's open
A unified description that handles both sectors explicitly at the Planck scale, not just as regime exits from each. Paper-6 Type IV exits (singularities, Planck discreteness, topology change) are acknowledged as loss of smooth structure, but there is no replacement for the smooth structure in that regime.

### Why this is not just "quantum gravity"
Standard quantum-gravity programs try to replace GR + QM with a deeper theory. APF's claim is that GR and QM are both regime-dependent representations of admissibility. The deeper layer is admissibility itself, not a new dynamical theory. So the "quantum-gravity" question inside APF is: what is the admissible class at the Planck scale? Not: what is the correct Planck-scale dynamics?

---

## 10. Matter-content bridge tightening (Paper 5 § 8)

### What's stated
Quantum structure and matter content are jointly constrained by the same capacity architecture (Paper 5, v2.0-PLEC, §8.1).

### What's open
The bridge from quantum-bookkeeping structure to specific gauge and matter content is presently a structural coexistence claim under a common capacity architecture. Tightening it to a direct single-theorem derivation ("Paper 5's admissibility forces Paper 4's 45 fermions") is open work. Paper 5 tags this as "mixed" in its epistemic-tag table.

### What would close it
A theorem chain from L_ST_algebra / L_ST_Hilbert (the Connes spectral triple machinery) to T_field (the 45-fermion count) without going through separate gauge-chain inputs.

---

## Meta: how to add a new entry to this list

If, in the course of working on APF, you discover a question that is (a) load-bearing for future work, (b) not currently resolved, and (c) not already catalogued here, update the canonical source at `__APF Library/Codebase/APF_Codebase_vX.X/OPEN_PROBLEMS.md` (or the relevant paper's Appendix C) and flag for inclusion in this template on the next rebuild of release repos.

Do not add entries to this file in a specific release — it will be overwritten. The canonical source is in the template folder: `__APF Library/skills/create-repo/templates/OPEN_PROBLEMS.md.template`.

---

*Generated by `create-repo`. Canonical at time of release. For the living workplan, see `__APF Library/APF Reference Docs/Reference - APF Paper Update Work Plan v2.md`.*
