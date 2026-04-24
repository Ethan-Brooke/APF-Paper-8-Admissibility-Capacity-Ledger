# Claims Ledger — Paper 8

Every load-bearing claim in Paper 8 maps here to a status tag, a proof location, a code check, and the failure mode a skeptical reviewer would attack. This file is the Paper-8 pilot for Phase 21.2 of the canonical work plan (AI-Onboarding Audit). If a claim appears in the supplement or main paper without a row here, that is a gap to be closed.

**Reading order.** Start with the rows tagged **nontrivial** — those are the structural core. Rows tagged **definitional**, **standard**, or **arithmetic** are consistency machinery around the core; they are not where the paper's weight sits.

**Status tag vocabulary** (sharper than `[P] / [C] / [P+lattice]`):

| Tag | Meaning |
|---|---|
| `nontrivial` | Load-bearing structural claim. Paper 8 stands or falls here. |
| `modest` | Definitional or near-definitional — true by construction given the record. |
| `standard` | Standard result from the mathematical literature; imported, not re-proved. |
| `arithmetic` | A numerical identity verifiable by direct calculation once the formula is chosen. |
| `empirical` | Comparison to observational data; not a derivation. |
| `conjectural` | Marked `[C]` — open, flagged in the paper, not claimed as proved. |
| `provenance` | The check reads a value from upstream and verifies it agrees with CANON. Cross-checks drift; does not derive the value. |

## Core claims

| # | Claim | Status | Local or imported? | Proof location | Code check | Failure mode |
|---|---|---|---|---|---|---|
| 1 | The ACC record is a two-scalar object $(K, \mathrm{ACC})$ with $\mathrm{ACC} = K \ln d_{\mathrm{eff}}$. | modest | local | Supp §2 | `ACC.__init__`, `ACC.value` in `apf/unification.py` | nonsensical record if $d_{\mathrm{eff}} \leq 0$; refused at construction. |
| 2 | Six regime projections $\pi_T, \pi_G, \pi_Q, \pi_F, \pi_C, \pi_A$ are well-defined on the ACC record. | modest | local | Supp §2 | `pi_T`, `pi_G`, `pi_Q`, `pi_F`, `pi_C`, `pi_A` in `apf/unification.py` | undefined projection raises; tested for each. |
| 3 | $I_1$ horizon-convention identity. | modest | local | Supp §D.1 | `_identity_I1_holographic`, `check_L_horizon_convention_equivalence` | wrong area convention; factor-of-4 error. |
| 4 | $I_2$ gauge-cosmological bridge. **The only nontrivial bridge claim of Paper 8.** | **nontrivial** | local formal kernel (Theorem 1.1) + imported irrep motivation (Paper 2 / Paper 4) | Supp §1 Theorem 1.1 + §D.2 integer-level corollary | `_identity_I2_gauge_cosmological`, `check_T_ACC_unification` | $V_\Lambda$ non-unique; $V_{\mathrm{local}}$ defined by what it excludes (circular); $V_{\mathrm{local}}$ irrep content unspecified. |
| 5 | $I_3$ thermo–quantum closure $S_{\mathrm{vN}}(\rho_{\max}) = \ln \dim \mathcal{H} = \mathrm{ACC}$. | standard | local proof, standard finite-dim result | Supp §D.3 | `_identity_I3_thermo_quantum` | breaks on non-max-mixed state; tied to finite-dim Hilbert assumption. |
| 6 | $I_4$ action–thermo closure $\lim_{\beta\to 0} \ln Z(\beta) = K \ln d_{\mathrm{eff}}$. | standard | local proof, standard partition-function argument | Supp §D.4 | `_identity_I4_action_thermo` | non-convergent $Z$ if interface not finite-dim; tested at $\beta = 0$ limit. |
| 7 | $I_2$ residual partition $3 + 16 + 42 = 61$. | provenance | imported from L_equip (cosmology.py) | Supp §1 + §E.3 | `acc_SM` consumes DAG values populated by `check_L_equip` | sector identification breaks (baryon = generations, dark includes Higgs); this is the attack surface. |
| 8 | Unique 42-dim $G_{\rm SM}$-invariant complement $V_\Lambda \subset V_{61}$. | nontrivial (paper proof complete; standalone witness landing Phase 22.2.a) | local formal kernel | Supp §1.4–1.5 Theorem 1.1 (Maschke); paper proof complete | `check_T_FormalKernel_VLambda_uniqueness` (landing Phase 22.2.a); downstream witness `check_T_interface_sector_bridge` in `gravity.py` | uniqueness depends on $V_{\rm local}$ being specified by irrep content, not by dimension alone — adversarial dimension-only-complement test in Phase 22.2.b. |
| 9 | $\Omega_\Lambda = \dim V_\Lambda / \dim V_{61} = 42/61$. | structural consequence | downstream of #4 + #7 | Supp §1.4 (iv) + §O.3 | `check_L_Omega_Lambda_is_entropy_fraction`, `pi_C` output | wrong if either partition (#7) or bridge (#4) fails. |
| 10 | $\Omega_b = 3/61, \Omega_c = 16/61$ as cosmological baryon + dark fractions. | provenance + arithmetic | downstream of #7 + matter-sector identification | Supp §1.3 + L_equip | `check_L_Omega_b_is_entropy_fraction`, `check_L_Omega_c_is_entropy_fraction` | fails if 3 or 16 is not the right sector count; this inherits the #7 attack surface. |
| 11 | $\rho_\Lambda / M_{\rm Pl}^4 = 42/102^{62}$ numerical formula. | arithmetic + conjectural (coefficient) | local formula composition | Supp §I + main §8 | `check_L_Lambda_absolute_numerical_formula` | exhaustive coefficient scan finds near-degenerate competitors (19/45 numerically closer); structural coefficient is `[C]`. |
| 12 | Numerical match $\rho_\Lambda / M_{\rm Pl}^4 \approx 10^{-122.91}$ agrees with Planck 2018 at 3% (= 8% in decades). | empirical | local arithmetic + observational input | Supp §I + main §8 | same as #11 | numerical agreement alone is insufficient; requires structural-coefficient derivation to be accepted. |
| 13 | $H_0 = 70.03$ km/s/Mpc derived from $\rho_\Lambda$ + GR $\rho_{\rm crit}$. | arithmetic | local algebraic inversion | Supp §I.3 | `check_T_Lambda_to_H0_inversion` | wrong GR form; H_0 tension framing is a falsifier, not a resolution. |
| 14 | $T_{\rm CMB} / M_{\rm Pl}$ parallel formula $(T_{\rm CMB}/M_{\rm Pl})^4 = 48/102^{64}$ matches FIRAS at 0.33%. | arithmetic + empirical | local formula | Supp §I + `thermal_absolute.py` | `check_T_T_CMB_absolute_formula` | same coefficient-degeneracy attack as #11. |
| 15 | Type III$_1$ factor at infinite volume. | standard | standard algebraic QSM (Araki 1968, Connes 1973, Kastler–Haag) | Supp §H.5 Theorem 7.14' | structural argument only; operator-level realisation tested in `quantum_operator_derivation.py` | assumptions A1 (non-degenerate local $H$), A2 (translation action), A3 (limit state) misstated. |
| 16 | SM ACC record uniqueness under four stated upstream identifications (U1)–(U4). | nontrivial (supplementary) | local | Supp §E.3 Theorem 5.6 | no standalone code check; uniqueness argument is prose | attacker replaces one of (U1)–(U4); theorem must exhibit contradiction. |
| 17 | Bayes-factor sanity check (prior-dependent 0.7 / 1.1 / 2.3 under flat / Gaussian-0.1 / Planck-matched priors). | empirical + arithmetic | local arithmetic | Supp §I'.4 | `check_L_Sigma_m_nu_suggestive`, `check_L_eta_does_not_fit_cleanly` (neighbouring checks) | **do not** read as full cosmological likelihood; it is a restricted sanity check. |
| 18 | Minimal working example at $K=3$, $d_{\rm eff}=4$: all six π projections verified by hand; one bridge identity verified; operator expectation verified. | modest | local pedagogical instance | Supp §R | `minimal_working_example/toy_interface_numpy.py` | toy mismatches would reveal a machinery bug; not relevant to SM-scale claims. |
| 19 | $C_{\mathrm{vacuum}} = 42$ derived from $27 + 3 + 12$ structural decomposition (T11). | nontrivial (independent route) | local, in `apf/gravity.py` T11 | Supp §T11 | `check_T11`, `check_L_self_exclusion` | are 27, 3, 12 independently specified, or do they trace to the same upstream as the matter-sector subtraction? This test is the second route to 42 and worth a direct cross-check. |

## Attack surface priority (for a hostile reviewer)

Claims tagged **nontrivial** and **provenance** are the priority attack surfaces. Rows 4, 7, 8, 10, 19 in that order. The coefficient-degeneracy concern on 11/12/14 is honestly disclosed as `[C]` in the paper; it is not a hidden vulnerability.

## How this ledger interacts with the anti-smuggling tests

`apf/test_no_smuggling.py` contains mutation tests that specifically probe claim 7 (residual partition provenance). As Phase 22 proceeds, expand coverage to claim 8 (Maschke uniqueness under $V_{\rm local}$ replacement), claim 11 (coefficient substitution 42/102 → 19/45 and confirm `[C]` status holds), and claim 19 (27 + 3 + 12 route independence from L_equip).

---

*Pilot rows: 19. Expected after full Phase 22 rollout: ~25–30. Anything in the supplement theorem index without a row here is an audit gap.*
