# Reviewers' Guide — APF Paper 8: The Admissibility-Capacity Ledger

This guide is the physics-first walkthrough for peer reviewers of The Admissibility-Capacity Ledger. It addresses the structural assumptions, anticipated objections, and falsification surfaces of the paper's argument. The intended audience is a peer reviewer with a background in physics or mathematics.

The guide complements three other verification routes:
- The [executable codebase](apf/) (run `python run_checks.py` for the full 36-theorem verification);
- The [Colab notebook](APF_Reviewer_Walkthrough.ipynb) (zero-install, theorem-by-theorem cells with prose + code);
- The [interactive derivation DAG](https://ethan-brooke.github.io/Admissibility-Capacity-Ledger/) (hover-for-details, click-for-dependencies, animated verification).

---

## The derivation in 36 steps

The full chain from the foundational commitments to this paper's results:

**Step 1 — ACC** (`ACC`)

  


**Step 2 — pi_T** (`pi_T`)

  Thermodynamic projection: S(rho, Gamma) = ACC(rho, Gamma).


**Step 3 — pi_G** (`pi_G`)

  Gravitational projection: S_BH = A/(4 ell_P^2) = ACC_horizon.


**Step 4 — pi_Q** (`pi_Q`)

  Quantum projection: dim H(rho, Gamma) = exp(ACC) = N.


**Step 5 — pi_F** (`pi_F`)

  Gauge/field projection: K(rho, Gamma) — the structural capacity.


**Step 6 — pi_C** (`pi_C`)

  Cosmological projection: residuals as fractions over K.


**Step 7 — pi_A** (`pi_A`)

  Action projection: Z(beta) = sum_g exp(-beta H(g)).


**Step 8 — _identity_I1_holographic** (`_identity_I1_holographic`)

  Identity I1 (Holographic): S_BH = ln(dim H_horizon) = ACC_horizon.


**Step 9 — _identity_I2_gauge_cosmological** (`_identity_I2_gauge_cosmological`)

  Identity I2 (Gauge-cosmological): K in pi_F equals denominator in pi_C.


**Step 10 — _identity_I3_thermo_quantum** (`_identity_I3_thermo_quantum`)

  Identity I3 (Thermo-quantum): S_vN(rho_max_mixed) = ln(dim H) = ACC.


**Step 11 — _identity_I4_action_thermo** (`_identity_I4_action_thermo`)

  Identity I4 (Action-thermo): lim_{beta -> 0} ln Z(beta) = ACC.


**Step 12 — T_ACC_unification** (`check_T_ACC_unification`)

  T_ACC_unification: Admissibility-Capacity Ledger unification [P].


**Step 13 — T_three_level_unification** (`check_T_three_level_unification`)

  T_three_level_unification — Three-level refinement of T_ACC [P]. (tier 4)


**Step 14 — T_subspace_functors_unified** (`check_T_subspace_functors_unified`)

  T_subspace_functors_unified [P] — All three 14c functors delivered. (tier 4)


**Step 15 — T_fractional_reading_equivalence** (`check_T_fractional_reading_equivalence`)

  T_fractional_reading_equivalence [P] — FRE at any ACC interface.


**Step 16 — L_Omega_Lambda_is_entropy_fraction** (`check_L_Omega_Lambda_is_entropy_fraction`)

  L_Omega_Lambda_is_entropy_fraction [P] — the dark-energy eureka.


**Step 17 — L_Omega_b_is_entropy_fraction** (`check_L_Omega_b_is_entropy_fraction`)

  L_Omega_b_is_entropy_fraction [P] — Baryon fraction as entropy allocation.


**Step 18 — L_Omega_c_is_entropy_fraction** (`check_L_Omega_c_is_entropy_fraction`)

  L_Omega_c_is_entropy_fraction [P] — CDM fraction as entropy allocation.


**Step 19 — T_residual_entropy_closure** (`check_T_residual_entropy_closure`)

  T_residual_entropy_closure [P] — Entropy fractions sum to 1.


**Step 20 — T_FRE_SM_to_entropy_dictionary** (`check_T_FRE_SM_to_entropy_dictionary`)

  T_FRE_SM_to_entropy_dictionary [P] — Full SM entropy dictionary. (tier 4)


**Step 21 — _identity_composition_tests** (`_identity_composition_tests`)

  


**Step 22 — L_Lambda_absolute_numerical_formula** (`check_L_Lambda_absolute_numerical_formula`)

  L_Lambda_absolute_numerical_formula [P] — APF coefficient formula matches Λ.


**Step 23 — T_Lambda_coefficient_degeneracy_audit** (`check_T_Lambda_coefficient_degeneracy_audit`)

  T_Lambda_coefficient_degeneracy_audit [C] — Numerical vs structural privilege.


**Step 24 — T_Lambda_absolute_extended_formula** (`check_T_Lambda_absolute_extended_formula`)

  T_Lambda_absolute_extended_formula [P] — Universal cosmological density formula.


**Step 25 — T_Lambda_to_H0_inversion** (`check_T_Lambda_to_H0_inversion`)

  T_Lambda_to_H0_inversion [P] — APF prediction for the Hubble constant.


**Step 26 — T_T_CMB_absolute_formula** (`check_T_T_CMB_absolute_formula`)

  T_T_CMB_absolute_formula [P] — APF predicts T_CMB via thermal C/d_eff^k.


**Step 27 — T_Lambda_absolute_bulletproof** (`check_T_Lambda_absolute_bulletproof`)

  T_Lambda_absolute_bulletproof [P] — Full robustness of Λ-absolute claim.


**Step 28 — _build_model** (`_build_model`)

  Build the (H_micro, H_total, P_vac_slot0) tuple for a model interface.


**Step 29 — T_Lambda_partition_function_at_beta_zero** (`check_T_Lambda_partition_function_at_beta_zero`)

  T_Lambda_partition_function_at_beta_zero [P] — Operator-level β → 0 limit.


**Step 30 — T_Lambda_d2_operator_derivation** (`check_T_Lambda_d2_operator_derivation`)

  T_Lambda_d2_operator_derivation [P over [P]+[P_structural]+[C]] —


**Step 31 — T_I3_operator_self_identification** (`check_T_I3_operator_self_identification`)

  T_I3_operator_self_identification [P_comp] — Composed top for I3.


**Step 32 — T_I4_operator_self_identification** (`check_T_I4_operator_self_identification`)

  T_I4_operator_self_identification [P_comp] — Composed top for I4.


**Step 33 — L_Sigma_m_nu_suggestive** (`check_L_Sigma_m_nu_suggestive`)

  L_Sigma_m_nu_suggestive [C, open] — neutrino mass sum consistent with k=15.


**Step 34 — L_eta_does_not_fit_cleanly** (`check_L_eta_does_not_fit_cleanly`)

  L_eta_does_not_fit_cleanly [C, open] — baryon-to-photon ratio scope.


**Step 35 — T_bridge_observer_independence_open** (`check_T_bridge_observer_independence_open`)

  T_bridge_observer_independence_open [C, open] —


**Step 36 — L_horizon_convention_equivalence** (`check_L_horizon_convention_equivalence`)

  L_horizon_convention_equivalence [P] — Conventions A and B agree.



Each step is a single theorem. The dependencies form a directed acyclic graph (DAG) with no cycles; the [interactive DAG](https://ethan-brooke.github.io/Admissibility-Capacity-Ledger/) shows multiple derivation paths to the same result, which is structural redundancy: many results are over-determined and would survive even if individual derivation steps were weakened.

---

## Pre-empting the operational-witness objection

A common objection to executable mathematical appendices is that the test-case witnesses (specific small examples used to verify the theorem) "smuggle in" the very structure being derived. Each test-case witness in this paper is the smallest concrete instance of a structure already proved in the abstract argument. Witnesses are constructed *after* the abstract claim is established, never before. The reviewer can re-derive each witness from the abstract derivation chain, not the other way around. If a witness appears to smuggle structure in, the corresponding `check_*` function will fail when the witness is replaced with an explicit countermodel.

---

## The structural assumptions

Beyond Axiom A1, this paper relies on the four PLEC components and a small set of structural regularity conditions:

- **A1 (Finite Enforcement Capacity)** — at any interface, the total enforcement cost across maintained distinctions is bounded.
- **MD (Minimum Distinction)** — every physical distinction has strictly positive enforcement cost (positive floor μ\* > 0).
- **A2 (Argmin selection)** — the realized configuration is the argmin of cost over the admissible set.
- **BW (Boundary Weight / Non-degeneracy)** — the cost spectrum is rich enough that the argmin is generically informative.

The four are pairwise logically independent (proved in Paper 1's Technical Supplement §1 with explicit countermodels). Together they define the variational structure named the **Principle of Least Enforcement Cost (PLEC)**: *reality is the minimum-cost expression of distinction compatible with finite capacity.*

_No paper-specific assumptions beyond the PLEC framework._

---

## The identifications

Several mappings between enforcement concepts and standard mathematical structures are used in the paper. They are **motivated, not deduced**, and reviewers should pay particular attention to their justification:

1. **Distinction ↔ binary partition of admissible configurations.** Spencer-Brown's distinction primitive (drawing the line between *this* and *not-this*) is identified with a binary measurement / partition at an enforcement interface. Motivated, not deduced.

2. **Cost ↔ resource expenditure.** The framework's $\varepsilon(d)$ functional is identified with the operational resource (energy / mass / barrier height) that maintains a distinction in any specific laboratory instance. This is a multi-resource correspondence; specific identifications appear in domain-specific papers (Paper 7 for action / Lagrangian, Paper 6 for spacetime metric, etc.).

3. **Capacity bound ↔ physical scarcity.** $C(\Gamma) < \infty$ at every interface is identified with the universally observed pattern that any specific physical system supports finitely many distinguishable states (Bekenstein--Hawking, Landauer, Ashby). Motivated by these traditions; not derived from them.

4. **Argmin ↔ realization.** The realized configuration at an interface is identified with the $\arg\min$ of total enforcement cost over the admissible set. This is a *locator*, not a process: nothing is doing the optimization; the argmin is what realization *means* under the framework.

5. **Hilbert space ↔ minimal representation of the noncommutative enforcement algebra.** Per the Paper 1 GNS argument, complex Hilbert space appears as the *unique* representational target for the algebraic structure that finite enforcement forces. Not assumed; selected via Frobenius.

These identifications are conventional rather than discovered: they are the choice that makes the framework's vocabulary and standard mathematics align. None is internally derivable.

---

## What is *not* proved in this paper

To prevent scope inflation, the following are flagged explicitly as *outside* the scope of Admissibility Capacity Ledger:

- **Specific gauge-group derivations** belong to Paper 2; this paper assumes them as imports if cited.
- **Specific particle content / generation count** belongs to Papers 2 and 4; not derived here.
- **Quantitative cosmological observables** beyond what is explicitly cited belong to Papers 3 (entropy / horizon) and 6 (geometry / DESI).
- **Quantum-gravity backreaction** is out of scope for any single paper in the current APF series; it is a future direction.
- **Numerical mass values** (absolute scales for $m_t$, $m_b$, $m_\nu$) are open problems noted in Paper 4 and the Engine paper; not within this paper's scope.
- **Spacetime dimension** (D = 4) is structurally derived in Paper 6 from Lovelock uniqueness; this paper assumes it where used.

If a reviewer concludes that the paper claims any of the above without supplying a proof for it, the reviewer is correct that the paper does not deliver that claim — those claims belong to other papers in the series and are explicitly flagged as such.

---

## How to falsify: attack surfaces

Each falsifier below targets a specific theorem or structural assumption. The corresponding code change in `apf/core.py` is also identified; reviewers can modify the codebase and re-run `python run_checks.py` to test each surface directly.

| # | Surface | Code-level test |
|---|---------|------------------|
| 1 | **Numerical disagreement.** A predicted observable disagrees with experiment beyond the framework's stated tolerance. | Modify the corresponding `check_*` to use the published experimental value; observe failure. |
| 2 | **Structural redundancy collapse.** Removing one PLEC component (A1, MD, A2, BW) leaves the derivation chain intact. | Comment out the test for the removed component in `apf/core.py`; observe other downstream checks failing. |
| 3 | **Reconstruction-program parity.** A standard quantum reconstruction (Hardy / CDP / Masanes--Müller) reaches the same conclusion with a strictly weaker assumption set. | Extract the GPT axioms used by the comparator; supply them as inputs to the relevant `apf/core.py` functions. |
| 4 | **Composition / locality break.** A multi-interface test with one interface deliberately violating $L_{\rm loc}$ does not produce the expected falsification mode. | Modify `check_L_loc` countermodel; rerun and observe expected vs. actual behavior. |
| 5 | **Cost-functional uniqueness fails.** An alternative cost functional satisfies all framework axioms equally well. | Replace the cost functional in `apf_utils.py`; observe whether downstream checks still pass. |
| 6 | **Scope-creep test.** A claim attributed to this paper is shown to actually require an unstated assumption. | Trace the claim's `\coderef`s through the bank dependency chain; identify any check that exits the paper's named scope. |

This is a structured threat model. If any of the surfaces fails empirically, the paper falsifies on that specific point.

---

## Reading the code

The codebase has three files in `apf/`:

- **`apf/core.py`** — the 36 theorem check functions for this paper. Each function constructs a mathematical witness, verifies the theorem's claim, and returns a structured result with name, dependencies, status, and key result.
- **`apf/apf_utils.py`** — exact arithmetic utilities (mostly `Fraction`-based; numpy/scipy where required by specific numerical lemmas).
- **`apf/bank.py`** — registry of all check functions in this repo, plus the `run_all()` runner.

Execution model: `run_checks.py` calls `bank.run_all()`, which iterates over every registered check and aggregates pass/fail/error counts. Individual checks can be invoked via `apf.bank.get_check('T_Born')()`.

---

## Scalar-to-matrix boundary

A characteristic feature of the APF program is that the early derivations use only finite sets and exact rational arithmetic (no matrices, no complex numbers, no linear algebra). Matrices first appear at the GNS construction (T2 in Paper 1), as the *minimal representation* of structure that the earlier scalar-only theorems prove must exist. This paper inherits the scalar-to-matrix transition from Paper 1's $T_2$. No new derivations of matrix structure occur in this paper; matrices appear where Paper 1 already established them.

This stratification is a deliberate methodological commitment: matrices are derived as representations of an already-proved abstract structure, not assumed as the substrate of the framework.

---

## The complex-field justification

The complex numbers as the unique admissible amplitude field is not a postulate but a derived selection (Paper 1 T2c, with proved exclusions $P_{\rm tom}$ and $P_{\rm cls}$ ruling out $\mathbb{R}$ and $\mathbb{H}$ respectively). This paper does not re-derive the complex-field selection; it inherits from Paper 1 ($T_{2c}$, with proved exclusions $P_{\rm tom}$ and $P_{\rm cls}$). For the original derivation, see Paper 1's [REVIEWERS_GUIDE.md](https://github.com/Ethan-Brooke/The-Enforceability-of-Distinction/blob/main/REVIEWERS_GUIDE.md).

---

## Citation and Zenodo

This repository is the executable mathematical appendix to APF Paper 8. The canonical archival deposit is at [https://doi.org/10.5281/zenodo.PENDING](https://doi.org/10.5281/zenodo.PENDING) (DOI: 10.5281/zenodo.PENDING).

---

*Generated by the APF `create-repo` skill on 2026-04-23. Codebase snapshot: v6.9 (frozen 2026-04-20; 355 verify_all checks, 342 bank-registered theorems, 48 quantitative predictions).*
