# The Admissibility-Capacity Ledger

### Interactive Mathematical Appendix to Paper 8 of the Admissibility Physics Framework

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19721384.svg)](https://doi.org/10.5281/zenodo.19721384) [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Ethan-Brooke/APF-Paper-8-Admissibility-Capacity-Ledger/blob/main/APF_Reviewer_Walkthrough.ipynb)

[Interactive Derivation DAG](https://ethan-brooke.github.io/APF-Paper-8-Admissibility-Capacity-Ledger/) · [Theorem Map](#theorem-mapping-table) · [Reviewers' Guide](REVIEWERS_GUIDE.md) · [The full APF corpus](#the-full-apf-corpus) · [Citation](#citation)

> **AI agents:** start with [`START_HERE.md`](START_HERE.md) — operational checklist that loads the framework context in 5–10 minutes. The corpus inventory and full file map are in [`ai_context/repo_map.json`](ai_context/repo_map.json).

---

## Why this codebase exists

The Admissibility-Capacity Ledger: a unifying record (K, d_eff) read through six regime projections (pi_T, pi_G, pi_Q, pi_F, pi_C, pi_A) with four consistency identities I_1-I_4. Paper 8's only non-trivial structural claim is the gauge-cosmological bridge I_2 (Theorem 1.1 of the Technical Supplement Formal Kernel): a unique G_SM-invariant 42-dimensional subspace V_Lambda of V_61 exists under the Standard-Model gauge action and T12 partition constraint, inducing the residual partition 3+16+42=61 and cosmological fraction Omega_Lambda = 42/61. I_1 is definitional under the cell-count convention K_horizon = A/(4 l_P^2); I_3 and I_4 are standard finite-dimensional closures. Supplement contains self-contained proof via Maschke semisimplicity + minimal working example at toy interface K=3 d_eff=4 auditable in ≤10 lines numpy + Planck 2018 multivariate Bayes-factor restricted sanity check + finite-volume C*-algebra formulation with Type III_1 classification at infinite volume. Four Phase 16-20 reviewer-response passes culminated in v2.9 main + v2.5 supplement.

This repository is the executable companion / verification layer for Paper 8's local claims and imported dependency subset.

The codebase is a faithful subset of the canonical APF codebase v6.9 (2026-04-23; 437 verify_all checks, 420 bank-registered theorems across 34 modules + `apf/standalone/`). Each manuscript result traces to a named `check_*` function in the bundled `apf/` subset — see the theorem-mapping table for module paths.

The codebase requires Python 3.8+ and NumPy / SciPy (some numerical lemmas use them; see `pyproject.toml`).

## How to verify

Three paths, in order of increasing friction:

**1. Colab notebook — zero install.** [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Ethan-Brooke/APF-Paper-8-Admissibility-Capacity-Ledger/blob/main/APF_Reviewer_Walkthrough.ipynb) Every key theorem is derived inline, with annotated cells you can inspect and modify. Run all cells — the full verification takes under a minute.

**2. Browser — zero install.** Open the [Interactive Derivation DAG](https://ethan-brooke.github.io/APF-Paper-8-Admissibility-Capacity-Ledger/). Explore the dependency graph. Hover any node for its mathematical statement, key result, and shortest derivation chain to A1. Click **Run Checks** to watch all theorems verify in topological order.

**3. Local execution.**

```bash
git clone https://github.com/Ethan-Brooke/APF-Paper-8-Admissibility-Capacity-Ledger.git
cd APF-Paper-8-Admissibility-Capacity-Ledger
pip install -e .
python run_checks.py
```

Expected output:

```
      Paper 8 (The Admissibility-Capacity Ledger): 36 passed, 0 failed, 36 total — verified in <minutes>
```

**4. Individual inspection.**

```python
from apf.bank import get_check
r = get_check('ACC')()
print(r['key_result'])
```

For reviewers, a [dedicated guide](REVIEWERS_GUIDE.md) walks through the logical architecture, the structural assumptions, and the anticipated objections.

---

## Theorem mapping table

This table maps every result in the manuscript to its executable verification.

| Check | Type | Summary |
|-------|------|---------|
| `ACC` | Other |  |
| `pi_T` | Other | Thermodynamic projection: S(rho, Gamma) = ACC(rho, Gamma). |
| `pi_G` | Other | Gravitational projection: S_BH = A/(4 ell_P^2) = ACC_horizon. |
| `pi_Q` | Other | Quantum projection: dim H(rho, Gamma) = exp(ACC) = N. |
| `pi_F` | Other | Gauge/field projection: K(rho, Gamma) — the structural capacity. |
| `pi_C` | Other | Cosmological projection: residuals as fractions over K. |
| `pi_A` | Other | Action projection: Z(beta) = sum_g exp(-beta H(g)). |
| `_identity_I1_holographic` | Other | Identity I1 (Holographic): S_BH = ln(dim H_horizon) = ACC_horizon. |
| `_identity_I2_gauge_cosmological` | Other | Identity I2 (Gauge-cosmological): K in pi_F equals denominator in pi_C. |
| `_identity_I3_thermo_quantum` | Other | Identity I3 (Thermo-quantum): S_vN(rho_max_mixed) = ln(dim H) = ACC. |
| `_identity_I4_action_thermo` | Other | Identity I4 (Action-thermo): lim_{beta -> 0} ln Z(beta) = ACC. |
| `check_T_ACC_unification` | Theorem | T_ACC_unification: Admissibility-Capacity Ledger unification [P]. |
| `check_T_three_level_unification` | Theorem | T_three_level_unification — Three-level refinement of T_ACC [P]. (tier 4) |
| `check_T_subspace_functors_unified` | Theorem | T_subspace_functors_unified [P] — All three 14c functors delivered. (tier 4) |
| `check_T_fractional_reading_equivalence` | Theorem | T_fractional_reading_equivalence [P] — FRE at any ACC interface. |
| `check_L_Omega_Lambda_is_entropy_fraction` | Lemma | L_Omega_Lambda_is_entropy_fraction [P] — the dark-energy eureka. |
| `check_L_Omega_b_is_entropy_fraction` | Lemma | L_Omega_b_is_entropy_fraction [P] — Baryon fraction as entropy allocation. |
| `check_L_Omega_c_is_entropy_fraction` | Lemma | L_Omega_c_is_entropy_fraction [P] — CDM fraction as entropy allocation. |
| `check_T_residual_entropy_closure` | Theorem | T_residual_entropy_closure [P] — Entropy fractions sum to 1. |
| `check_T_FRE_SM_to_entropy_dictionary` | Theorem | T_FRE_SM_to_entropy_dictionary [P] — Full SM entropy dictionary. (tier 4) |
| `_identity_composition_tests` | Other |  |
| `check_L_Lambda_absolute_numerical_formula` | Lemma | L_Lambda_absolute_numerical_formula [P] — APF coefficient formula matches Λ. |
| `check_T_Lambda_coefficient_degeneracy_audit` | Theorem | T_Lambda_coefficient_degeneracy_audit [C] — Numerical vs structural privilege. |
| `check_T_Lambda_absolute_extended_formula` | Theorem | T_Lambda_absolute_extended_formula [P] — Universal cosmological density formula. |
| `check_T_Lambda_to_H0_inversion` | Theorem | T_Lambda_to_H0_inversion [P] — APF prediction for the Hubble constant. |
| `check_T_T_CMB_absolute_formula` | Theorem | T_T_CMB_absolute_formula [P] — APF predicts T_CMB via thermal C/d_eff^k. |
| `check_T_Lambda_absolute_bulletproof` | Theorem | T_Lambda_absolute_bulletproof [P] — Full robustness of Λ-absolute claim. |
| `_build_model` | Other | Build the (H_micro, H_total, P_vac_slot0) tuple for a model interface. |
| `check_T_Lambda_partition_function_at_beta_zero` | Theorem | T_Lambda_partition_function_at_beta_zero [P] — Operator-level β → 0 limit. |
| `check_T_Lambda_d2_operator_derivation` | Theorem | T_Lambda_d2_operator_derivation [P over [P]+[P_structural]+[C]] — |
| `check_T_I3_operator_self_identification` | Theorem | T_I3_operator_self_identification [P_comp] — Composed top for I3. |
| `check_T_I4_operator_self_identification` | Theorem | T_I4_operator_self_identification [P_comp] — Composed top for I4. |
| `check_L_Sigma_m_nu_suggestive` | Lemma | L_Sigma_m_nu_suggestive [C, open] — neutrino mass sum consistent with k=15. |
| `check_L_eta_does_not_fit_cleanly` | Lemma | L_eta_does_not_fit_cleanly [C, open] — baryon-to-photon ratio scope. |
| `check_T_bridge_observer_independence_open` | Theorem | T_bridge_observer_independence_open [C, open] — |
| `check_L_horizon_convention_equivalence` | Lemma | L_horizon_convention_equivalence [P] — Conventions A and B agree. |

Each manuscript claim traces to a `check_*` function in the bundled `apf/` subset; modules include `core.py`, `gauge.py`, `gravity.py`, `unification.py`, `unification_three_levels.py`, `subspace_functors.py`, `fractional_reading.py`, `lambda_absolute.py`, `lambda_operator_derivation.py`, and others. See `ai_context/theorems.json` for the full module-to-check mapping. Every function can be called independently and returns a structured result including dependencies and verified content.

---

## The derivation chain

```
  Level 0: ACC · pi_T · pi_G · pi_Q · pi_F · pi_C · pi_A · _identity_I1_holographic · _identity_I2_gauge_cosmological · _identity_I3_thermo_quantum · _identity_I4_action_thermo · T_ACC_unification · T_three_level_unification · T_subspace_functors_unified · T_fractional_reading_equivalence · L_Omega_Lambda_is_entropy_fraction · L_Omega_b_is_entropy_fraction · L_Omega_c_is_entropy_fraction · T_residual_entropy_closure · _identity_composition_tests · L_Lambda_absolute_numerical_formula · T_Lambda_coefficient_degeneracy_audit · T_Lambda_absolute_extended_formula · T_Lambda_to_H0_inversion · T_T_CMB_absolute_formula · _build_model · T_Lambda_partition_function_at_beta_zero · T_I3_operator_self_identification · L_Sigma_m_nu_suggestive · L_eta_does_not_fit_cleanly · T_bridge_observer_independence_open · L_horizon_convention_equivalence
  Level 1: T_FRE_SM_to_entropy_dictionary · T_Lambda_absolute_bulletproof · T_Lambda_d2_operator_derivation · T_I4_operator_self_identification
```

The [interactive DAG](https://ethan-brooke.github.io/APF-Paper-8-Admissibility-Capacity-Ledger/) shows the full graph with hover details and animated verification.

---

## Repository structure

```
├── README.md                              ← you are here
├── START_HERE.md                          ← AI operational checklist; read-first for AI agents
├── REVIEWERS_GUIDE.md                     ← physics-first walkthrough for peer reviewers
├── interactive_dag.html                   ← interactive D3.js derivation DAG (also served at docs/ via GitHub Pages)
├── repo_map.json                          ← machine-readable map of this repo (root copy of ai_context/repo_map.json)
├── theorems.json                          ← theorem catalog (root copy of ai_context/theorems.json)
├── derivation_graph.json                  ← theorem DAG as JSON (root copy of ai_context/derivation_graph.json)
├── ai_context/                            ← AI onboarding pack (corpus map, theorems, glossary, etc.)
│   ├── AGENTS.md                          ← authoritative entry point for AI agents
│   ├── FRAMEWORK_OVERVIEW.md              ← APF in 5 minutes
│   ├── GLOSSARY.md                        ← axioms, PLEC primitives, epistemic tags
│   ├── AUDIT_DISCIPLINE.md                ← engagement posture for critique/proposal
│   ├── OPEN_PROBLEMS.md                   ← catalog of open problems + verdicts
│   ├── repo_map.json                      ← machine-readable map of this repo
│   ├── theorems.json                      ← machine-readable theorem catalog
│   ├── derivation_graph.json              ← theorem DAG as JSON
│   └── wiki/                              ← bundled APF wiki (concepts, papers, codebase)
├── apf/
│   ├── core.py                            ← A1 + PLEC (foundational)
│   ├── apf_utils.py                       ← exact arithmetic + helpers
│   └── bank.py                            ← registry and runner
├── docs/
│   └── index.html                         ← interactive derivation DAG (GitHub Pages)
├── APF_Reviewer_Walkthrough.ipynb         ← Colab notebook
├── run_checks.py                          ← convenience entry point
├── pyproject.toml                         ← package metadata
├── zenodo.json                            ← archival metadata
├── Paper_8_Admissibility_Capacity_Ledger_v2.9.tex                ← the paper
├── Paper_8_Admissibility_Capacity_Ledger_Supplement_v2.5.tex                ← Technical Supplement

└── LICENSE                                ← MIT
```

---

## What this paper derives and what it does not

**Derived:** (see Theorem mapping table above)

**Not derived here:** Specific results outside this paper's scope live in companion papers — see the corpus table below for the full 9-paper series.

---

## Citation

```bibtex
@software{apf-paper8,
  title   = {The Admissibility-Capacity Ledger},
  author  = {Brooke, Ethan},
  year    = {2026},
  doi     = {10.5281/zenodo.19721384},
  url     = {https://github.com/Ethan-Brooke/APF-Paper-8-Admissibility-Capacity-Ledger}
}
```

For the full citation lineage (concept-DOI vs version-DOI, related identifiers, bibtex for all corpus papers), see [`ai_context/CITING.md`](ai_context/CITING.md).

---

## The full APF corpus

This repository is **one paper-companion** in a 9-paper series. Each paper has its own companion repo following this same layout. The full corpus, with canonical references:

| # | Title | Zenodo DOI | GitHub repo | Status |
|---|---|---|---|---|
| 0 | What Physics Permits | [10.5281/zenodo.18439523](https://doi.org/10.5281/zenodo.18439523) | [`APF-Paper-0-What-Physics-Permits`](https://github.com/Ethan-Brooke/APF-Paper-0-What-Physics-Permits) | public |
| 1 | The Enforceability of Distinction | [10.5281/zenodo.18439200](https://doi.org/10.5281/zenodo.18439200) | [`APF-Paper-1-The-Enforceability-of-Distinction`](https://github.com/Ethan-Brooke/APF-Paper-1-The-Enforceability-of-Distinction) | public |
| 2 | The Structure of Admissible Physics | [10.5281/zenodo.18439274](https://doi.org/10.5281/zenodo.18439274) | [`APF-Paper-2-The-Structure-of-Admissible-Physics`](https://github.com/Ethan-Brooke/APF-Paper-2-The-Structure-of-Admissible-Physics) | public |
| 3 | Ledgers | [10.5281/zenodo.18439363](https://doi.org/10.5281/zenodo.18439363) | [`APF-Paper-3-Ledgers-Entropy-Time-Cost`](https://github.com/Ethan-Brooke/APF-Paper-3-Ledgers-Entropy-Time-Cost) | public |
| 4 | Admissibility Constraints and Structural Saturation | [10.5281/zenodo.18439397](https://doi.org/10.5281/zenodo.18439397) | [`APF-Paper-4-Admissibility-Constraints-Field-Content`](https://github.com/Ethan-Brooke/APF-Paper-4-Admissibility-Constraints-Field-Content) | public |
| 5 | Quantum Structure from Finite Enforceability | [10.5281/zenodo.18439433](https://doi.org/10.5281/zenodo.18439433) | [`APF-Paper-5-Quantum-Structure-Hilbert-Born`](https://github.com/Ethan-Brooke/APF-Paper-5-Quantum-Structure-Hilbert-Born) | public |
| 6 | Dynamics and Geometry as Optimal Admissible Reallocation | [10.5281/zenodo.18439445](https://doi.org/10.5281/zenodo.18439445) | [`APF-Paper-6-Dynamics-Geometry-Spacetime-Gravity`](https://github.com/Ethan-Brooke/APF-Paper-6-Dynamics-Geometry-Spacetime-Gravity) | public |
| 7 | Action, Internalization, and the Lagrangian | [10.5281/zenodo.18439513](https://doi.org/10.5281/zenodo.18439513) | [`APF-Paper-7-Action-Internalization-Lagrangian`](https://github.com/Ethan-Brooke/APF-Paper-7-Action-Internalization-Lagrangian) | public |
| 13 | The Minimal Admissibility Core | [10.5281/zenodo.18361446](https://doi.org/10.5281/zenodo.18361446) | [`APF-Paper-13-The-Minimal-Admissibility-Core`](https://github.com/Ethan-Brooke/APF-Paper-13-The-Minimal-Admissibility-Core) | public |
| — | Canonical codebase (v6.9) | [10.5281/zenodo.18529115](https://doi.org/10.5281/zenodo.18529115) | [`APF-Codebase`](https://github.com/Ethan-Brooke/APF-Codebase) | pending |

The canonical computational engine — the full bank of 342 theorems across 19 modules — is the **APF Codebase** ([Zenodo](https://doi.org/10.5281/zenodo.18529115)). Every per-paper repo is a faithful subset of that engine.

---

## License

MIT. See [LICENSE](LICENSE).

---

*Generated by the APF `create-repo` skill on 2026-04-23. Codebase snapshot: v6.9 (frozen 2026-04-20; 355 verify_all checks, 342 bank-registered theorems, 48 quantitative predictions).*
