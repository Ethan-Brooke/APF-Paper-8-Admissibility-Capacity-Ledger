---
type: paper
domain: apf
layer: 7-action
created: 2026-04-14
updated: 2026-04-21 (Phase 12.4 I4 action-thermo ingest: v2.0/v1.0 → v2.1/v1.1)
sources: []
---

# Paper 7 - Action
## "Action Principle, Internalization, and BSM Physics"

**Status:** **v2.1 main (23pp) + Technical Supplement v1.1 (19pp)** (Phase 12.4 I4 action-thermo ingest, 2026-04-21). Paper 7 is the receiver of T_ACC identity $I4$ (action-thermo: $\lim_{\beta\to 0}\ln Z = \mathrm{ACC}$) since Paper 7's Part II ($Z_0$ partition function) and Part III (spectral action = partition function via $L_{\mathrm{spectral\_action\_internal}}$) are precisely what $I4$ abstracts into the T_ACC ledger at the action-regime projection $\pi_A$. Earlier: v2.0/v1.0 (2026-04-17 full-scope rebuild).

**Changelog:**
- **v2.1 main + Supplement v1.1 (2026-04-21):** Phase 12.4 I4 action-thermo ingest. Supplement: new §sec:I4_register subsection at end of §sec:Z0 with Theorem "Paper 7 content of identity I4" proving $\lim_{\beta\to 0}\ln Z_0(\beta) = C_{\mathrm{total}}\ln(1+d_{\mathrm{eff}})$ and reconciling with $\mathrm{ACC} = K\ln d_{\mathrm{eff}}$ up to $O(\ln(1+1/d_{\mathrm{eff}}))$ correction; two follow-on remarks (Paper 7 as the $\pi_A$ regime factory for the T_ACC ledger; $I4$ as limit-form of the Ledger-Accounting Principle). Spectral-action pivot remark (§sec:SAeqZ) extended with T_ACC reading: spectral action = $\pi_A$ at the spectral-triple interface; $\Lambda^2 = 1/\beta$ correspondence pins the Boltzmann cutoff to the $I4$ trajectory; triple reading of $a_0 = 61 = C_{\mathrm{total}}$ (Seeley–DeWitt / $\pi_F$ / $\pi_C$). Main: Part II opening (§sec:partfn) T_ACC orientation paragraph identifying $Z(\beta)$ as the canonical $\pi_A$ realization; Part III opening (§sec:specaction) T_ACC reading paragraph identifying spectral-action = partition-function as the Boltzmann-cutoff specialization of $I4$. Compile-debug: `\end{remark}` doubled-close fixed; `\Ctot` macro defined only in supplement accidentally used in main (replaced with `C_{\mathrm{total}}`); stale cross-doc `\ref{thm:I4_P7}` replaced with prose pointer. Archive: `Old/Paper_7_*_pre-Phase12.4.{tex,pdf}`. Upstream reference: Three-Layer Ontology doc §7 (2ε-per-channel) and §8 (K and ACC two-scalar).
- **v2.0 live (2026-04-17)** — full-scope rebuild with Technical Supplement. **Main paper 22pp (was 21pp pre-wayfinder)**, supplement `Paper_7_Technical_Supplement_v1.0.pdf` (17pp), both compiled cleanly with pdflatex. Source of truth lives in the supplement; main paper is argument-first with ≤20-line proof sketches pointing to supplement theorems. **PLEC propagation and partition-function insertion are now native to the rebuild** — no longer source-pending (was Phase 2.8b under previous accounting). As of 2026-04-18, includes the series-wide "Operational primitives" wayfinder tcolorbox at start of Part I mapping five primitives (quantum of action ε*, partition function Z, saturation factorisation, spectral triple, Seeley-DeWitt expansion) to their definition and formalization sites.

**Prior version:** Feb-2026 6pp PDF `Paper 7 Minimal Quantum Action.pdf` (Zenodo 18604875) is preserved in `Old/Paper_7_Minimal_Quantum_Action_v1.0_Feb2026.pdf` and remains at the top level under its original filename.

**Live cross-references from other papers (satisfied by v2.0):** Paper 3 v3.0 §5 ledger-accounting principle now has its thermodynamic-formalism backing in Paper 7 §S.3 (saturation factorisation). Paper 3 v3.0 Appendix C "Does Not Claim" pointer to Paper 7 is also satisfied.

**Core claim:** Zero action is structurally inadmissible for record-forming paths; the minimum quantum of action is the MD component of PLEC (floor `ε*`). Under a spectral-triple internalization of the APF substrate, the Connes spectral action equals the APF partition function at Boltzmann cutoff, and the Seeley–DeWitt expansion of the resulting heat trace yields the Standard Model Lagrangian (Yang–Mills + Higgs) plus gravity (cosmological constant `a_0 = C_total = 61` + Einstein–Hilbert). BSM mechanisms (dark matter, inflation, baryogenesis) are named as structurally compatible without quantitative claims.

**Paper structure (five Parts):**
- **Part I — Quantum of Action** (§§2–6): primitives; paths/records with BFS; minimum-action Theorem (Isolation of Zero); multi-record `k!` factor
- **Part II — Partition Function** (§§7–9): `Z_0(β) = (1 + d_eff e^{-βε*})^{C_total}`; gapped regime `kT ≪ ε*` with 10^32 Planck margin; saturation factorisation `Z_sat = N_micro · exp(-β E_sat(η))`
- **Part III — Spectral Action = Partition Function** (§§10–12): Connes = APF identification; Seeley–DeWitt `a_0=61, a_2≈21.985, a_4≈87.201`; ½ coefficient from KO-dim 6; scalar-potential form with mass gap
- **Part IV — Geometric Internalization** (§§13–17): Kolmogorov continuum, chartability, Lovelock → Einstein unique in d=4, Coleman–Mandula, HKM causal order, Malament full-metric
- **Part V — SM Lagrangian Extensions** (§§18–22): one-loop β coefficients `b_1=41/10, b_2=-19/6, b_3=-7` from APF field content; Type-I seesaw; Pauli–Jordan spin-statistics; McKean–Singer index; Tannaka–Krein → SU(3)×SU(2)×U(1); ε ↔ ℏ identification table

Plus §23 falsifiers F1–F7, §24 summary with implications table, §25 conceptual status. Appendices: A (cost-parameter crosswalk), B (toy model), C (PLEC assumption inventory), D (changelog).

**Key theorems (all with `\coderef`):**
- `L_spectral_action_internal` (`check_spectral_action_internal`, `internalization.py`) — **structural pivot**: Connes spectral action = APF partition function at Boltzmann cutoff
- `L_normalization_coefficient` (`internalization.py`) — ½ on Tr(κ†κ) from KO-dim 6
- `L_scalar_potential_form` (`internalization.py`) — V(H,σ) form from spectral invariance + A1 + K3
- `L_kolmogorov_internal`, `L_chartability`, `L_lovelock_internal`, `L_coleman_mandula_internal` (`internalization_geo.py`)
- `L_HKM_causal_geometry`, `L_Malament_uniqueness` (`extensions.py`) — internalized here rather than forward-referencing Paper 6 (stand-alone convention)
- `T6B_beta_one_loop`, `L_seesaw_type_I`, `L_Pauli_Jordan`, `L_McKean_Singer_internal`, `L_Tannaka_Krein` (`extensions.py`)
- Minimum-action theorems (Theorem 5.2 Isolation of Zero; Theorem 6.1 `k!` multi-record factor) — proved stand-alone in supplement §S.1

**Supplement structure (Paper 1/2 mirror):** front-matter PLEC box + what-is-assumed table + framing convention. §S.1 minimum action (4 theorems). §S.2 partition function (3 theorems). §S.3 saturation factorisation (Theorem 3.2 + 2 hedging remarks). §S.4 spectral action = partition function (Theorem 4.2 equivalence, 4.3 Seeley–DeWitt, 4.4 uniqueness). §S.5 geometric internalization (6 theorems with hypotheses H2 Lipschitz and H3 HKM flagged). §S.6 SM Lagrangian extensions (5 theorems). §S.7 red-team (3 additional hypotheses audited; wiki parked-material hedges adopted explicitly). Appendix A theorem index with code references and epistemic tags. Appendix B dependency diagram.

**Hypotheses flagged in red-team §S.7:**
- **H1 Spectral-triple ansatz** — Assumption 4.1: APF algebra + Hilbert carry self-adjoint Dirac `D` with KO-dim 6. Additional constitutive structure; not derivable from PLEC alone.
- **H2 Lipschitz cost regularity** — used for Kolmogorov continuum + chartability. Analogous to MD in character.
- **H3 HKM future-distinguishing / density** — standard Lorentzian-causality hypothesis; used for causal-order uniqueness and metric recovery.

Each hypothesis is marked at point of use and the theorem index marks downstream results as `[P+H1/H2/H3]` rather than `[P]`.

**Hedges from wiki parked-material block (adopted in §S.7):**
- Fermionic-substrate reading of APF matter content: "suggestive, not conclusive"
- β- and η-independence at saturation: framed as clean factorisation, not as independent theorems
- "Λ·G is exact" contest line: replaced by "Λ·G forced by counting, not tuned"

**Implementation:** `internalization.py` + `internalization_geo.py` + `extensions.py` (14 referenced checks)

**Status notes:**
- v2.0 main + supplement shipped 2026-04-17
- PLEC propagation + partition-function insertion now native
- BSM predictions: mechanisms named without quantitative claims per Q2=i conservative scope

---

## Parked material: Partition-function / thermodynamic formulation (for future Paper 7 expansion)

*Added 2026-04-17 following review of Ethan's partition-function investigation. Decision: this material belongs in Paper 7, not Paper 1 or 2. Full review and scope rationale in [[Log#2026-04-17 review]].*

**What to include (with hedges):**

- **Non-interacting partition function.** `Z_0(β) = (1 + d_eff · e^{-β ε*})^{C_total}`, derived as a sum over subsets of the C_total capacity types with uniform gap ε*. At d_eff = 1, recovers Fermi-Dirac form for binary modes; at d_eff = 102, generalises to multi-level (Potts-like) occupation.
- **Gapped-ground-state framing.** APF is valid in the regime `kT ≪ ε*`. If ε* ~ E_P, this covers the observable universe with margin ~10^32 and breaks down at the Planck era — consistent with `L_singularity_resolution` and the inflation-staircase picture. This is Paper 7's useful framing result for APF's regime of validity.
- **Saturation factorisation.** `Z_sat(β, η) = N_micro · exp(-β E_sat(η))` with `N_micro = d_eff^{C_total}` and `E_sat(η) = C_total ε* + Σ_{i<j} η(i,j)`. States cleanly the structural/thermodynamic separation: combinatorial content lives in `N_micro`, thermodynamic content in the energy shift. One paragraph in Paper 7 supplement is enough; it answers "what about loop/thermal/interaction corrections to S_dS" in a line.

**What to NOT include (or include only with heavy hedging):**

- **"Fermionic substrate" interpretation of SM matter content.** The partition function has the same mathematical form as non-interacting fermions at d_eff = 1, but "capacity committed vs. not" is not the same ontological category as "Pauli-excluded single-particle state." The report (`APF_Partition_Function_Report.md`) hedges this correctly ("suggestive rather than conclusive") — keep that hedge. The stronger claim ("fermionic character of SM matter is not an accident") is overreach on the current argument.
- **"Saturation exactness" framed as a theorem.** β-independence and η-independence of S_dS at saturation are mathematically correct but near-trivial: once you condition on a single mode configuration, multiplicity is a count and counts don't depend on energy scales. `Σ η(i,j)` on a single configuration is a constant energy shift, which by construction doesn't change multiplicity. Do not publish these as independent named theorems. The entire content fits in one sentence: *"In the saturated ensemble, Z_sat factorises into a β- and η-independent multiplicity and a Boltzmann factor on the fixed saturation energy; therefore S_dS depends only on the counting structure."*
- **"Λ·G is exact, not approximate" contest/marketing line.** This invites the standard "tell me about your corrections" attack and misdirects from where the actual derivational weight lives (`L_count` → C_total = 61; `L_self_exclusion` → d_eff = 102; holographic relation → Λ·G = 3π · exp(-S_dS)). Safer equivalent: *"Λ·G is forced by counting, not tuned. Each factor in the count is itself structurally forced."*

**Open research direction (beyond Paper 7's current scope):**

- **Sub-saturated regime.** Saturation (k = C_total = 61) is the fixed endpoint where the partition function structure is trivial. Interactions η(i,j) and thermal β-dependence do real work in the sub-saturated regime (pre-inflationary trajectory k=1 → k=61). If APF has observable early-universe signatures beyond `L_singularity_resolution`, they live there. Post-contest research direction.
- **d_eff = 102 generalisation.** If each mode supports 102 internal states, the "fermionic" Fermi-Dirac form is replaced by a Potts-like (1 + 101 e^{-βε*})^61. Whether this gives a unified thermodynamic derivation of Bekenstein entropy, de Sitter entropy, and cosmological constant depends on what d_eff = 102 represents structurally. Lead worth following after contest.

**Supporting artifacts (user working notes, not checked into APF Library):**

- `APF_Partition_Function_Report.md` — original investigation of Z(β) and β's physical meaning
- `beta_physical_meaning.py` — numerical exploration of β vs. epoch k
- `APF_Saturation_Exactness_Report.md` — follow-up on β- and η-exactness claims
- `finite_temp.py` — thermal corrections calculation (Schottky peak at β ≈ 4.62)
- `interaction_corrections.py` — pairwise interaction corrections, 4-mode verification

## See also
- [[Derivation Chain]]
- [[Open Problems]]
- [[Paper 4 - Constraints]]
- [[Cosmological Predictions]]
- [[T_ACC Unification]] — Paper 7 supplies the $\pi_A$ projection; $I4$ identity realized here
- [[Paper 3 - Ledgers]] — ACC formalization home (Supplement v1.2 §7)
