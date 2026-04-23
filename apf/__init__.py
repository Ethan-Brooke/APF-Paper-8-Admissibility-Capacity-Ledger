"""APF v6.9 — Admissibility Physics Framework.

347 bank-registered theorems across 20 modules (plus apf/standalone/);
362 total verify_all checks. Zero postulates.

v6.9 (2026-04-20 late evening): interface-sector bridge pass.
  New check_T_interface_sector_bridge (tier 4, [P]) and auxiliary
  check_L_global_interface_is_horizon (tier 3, [P]) in apf/gravity.py.
  The T12 interface partition V_61 = V_local (+) V_global governs the
  T_horizon_reciprocity second-epsilon sector decomposition:
    |Sector A| = |V_61 \\ {self}| = 60,
    |Sector B| = dim V_global    = 42,
    d_eff      = 60 + 42         = 102.
  The auxiliary lemma identifies T12's non-finite-interface stratum with
  the dS horizon absorber subspace (T_Bek), closing Step 4 of the bridge
  theorem without an implicit cross-module identification. Promotes the
  two "42"s in T11 and L_self_exclusion to a single witnessed geometric
  identity on V_global. Strictly finer than T_ACC identity I2 (integer K
  → subspace V_global). Bank 345-loaded → 347 (matches EXPECTED_THEOREM_COUNT
  = 347; absorbs pre-existing silent -2 drift); verify_all 360 → 362.

v6.9 (2026-04-20 afternoon): T_ACC unification pass.
  New apf/unification.py module — the Admissibility-Capacity Ledger.
  ACC record (two scalars K, ACC = K ln d_eff), six regime projections
  (pi_T, pi_G, pi_Q, pi_F, pi_C, pi_A), three canonical-interface factories
  (acc_SM, acc_horizon, acc_quantum), five bank-registered checks:
    I1_holographic           — S_BH = ln(dim H_horizon) = ACC_horizon.
    I2_gauge_cosmological    — K (pi_F at SM) = K (pi_C denominator at SM).
    I3_thermo_quantum        — S_vN(rho_max) = ln dim H = ACC.
    I4_action_thermo         — ln Z(beta) -> ACC as beta -> 0.
    T_ACC_unification        — composed theorem re-running I1..I4.
  Bank 342 -> 347; verify_all 355 -> 360; modules 19 -> 20.

v6.9 (2026-04-18): PLEC formalization.
  New apf/plec.py module with Regime R + five-type regime-exit taxonomy:
    Regime_R                 — R1..R4 joint validity; PLEC well-posedness.
    Regime_exit_Type_I       — collapse of admissible variation (saturation).
    Regime_exit_Type_II      — minimizer nonuniqueness (branching).
    Regime_exit_Type_III     — change of admissible class (record locking).
    Regime_exit_Type_IV      — loss of smooth / local structure.
    Regime_exit_Type_V       — pure representational redundancy.
  New A9_closure in apf/gravity.py unifying the Lovelock prerequisites
    A9.1..A9.5 (locality, covariance, conservation, second-order,
    propagation) dispersed across core/gravity/spacetime/internalization_geo.
  Papers 5/6 v2.0-PLEC now code-anchored (coderef pass complete).

v6.8: Canonicalization (2026-04-18).
  335 bank-registered / 348 verify_all; 18 modules; archive + naming discipline.

v6.7: Option 3 Work Plan — Phases 1–6 complete.
  Phase 1 (seesaw gap): L_seesaw_from_A1 [P] — 9-link chain, zero imports.
  Phase 2 (mass matrix): L_mass_from_capacity [P] — 11-link chain, zero FN.
    L_multiplicative_amplitude, L_Yukawa_bilinear, RT_FN_vs_capacity.
  Phase 3 (texture): L_texture_from_capacity [P] — 10-link chain, zero Fritzsch.
    L_GJ_from_capacity, RT_texture_chain.
  Phase 4 (bridges): L_bridges_closed [P] — all 5 bridges now theorems.
    RT_bridge_audit.
  Phase 5 (Theorem R): RT_R1/R2/R3 [P] — not circular, R3 rewritten.
  Phase 6 (NCG): L_NCG_status [P] — 11/11 items derived, zero physics imports.
    RT_NCG_no_physics_import. Long-term: derive formalism itself (math research).

v6.6: δ_PMNS resolution + DESI DR2 cosmological confrontation.
  L_seesaw_factorization, L_PMNS_CP_corrected, L_DESI_DR2_confrontation,
  L_joint_cosmo_neutrino, L_top_mass_hint.

v6.5: Down-sector closure + NNLO precision.
  L_Higgs_curvature_channel, L_NNLO_Fritzsch, L_sin2_oneloop, L_lepton_GJ.

v6.4: Gauge uniqueness + Weinberg angle Cauchy + CKM resolution.
"""
__version__ = '6.9'
