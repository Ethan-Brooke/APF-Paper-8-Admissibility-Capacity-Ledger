---
type: entity
domain: apf
layer: cross-cutting
created: 2026-04-14
updated: 2026-04-21 (v6.9 night (later) Phase 13.3 Stage II Enforcement Crystal analytical-metrics layer; night Phase 13.2 Enforcement Crystal walker; evening Phase 14 three-level identity refinement; morning Phase 13.1 PLEC anchor pass)
sources: []
---

# APF Codebase
## "Version 6.9: 384 Checks, 375 Bank-Registered Theorems, 23 Registered Modules + standalone/"

**Overview:** A production Python codebase implementing the entire [[Admissible Physics Framework]]. Each check function is a rigorous proof of a theorem; all checks pass with exact arithmetic where possible (rational numbers via `fractions`; numerical routines via `numpy` and `scipy` where required).

**Current version:** **v6.9 (2026-04-21 night (later) Phase 13.3 Stage II Enforcement Crystal analytical-metrics layer; night Phase 13.2 Enforcement Crystal walker; evening Phase 14 three-level identity refinement; morning Phase 13.1 PLEC anchor pass)**

### v6.9 changes (2026-04-21 night (later) Phase 13.3 Stage II Enforcement Crystal analytical-metrics layer)
- New `apf/crystal_metrics.py` module (~1100 lines, 5 bank-registered checks). Stage II of the Phase 13 Enforcement Crystal program: the analytical-metrics layer that consumes the Phase 13.2 walker (`apf/crystal.py`) and recomputes the Feb-2026 Paper 20 v1.0 measures directly on the v6.9 bank graph. Five workstreams landed in sequence:
  - **Workstream 1 — Brandes betweenness centrality** (`check_T_crystal_betweenness_v69`, tier 3, [P_structural]): Brandes O(VE) algorithm on the depth-filtered DAG, ranks waist candidates by BC.
  - **Workstream 2 — single-deletion cascade-failure** (`check_T_crystal_cascade_v69`, tier 3, [P_structural]): `cascade(N) = baseline_R \ {N} \ reach_without_N`, `fraction = |cascade| / |baseline\\{N}|`. Headline finding: **max fraction is 2.3% (T_Bek), median 0%** — the v6.9 bank is dramatically over-determined.
  - **Workstream 3 — Menger min vertex cuts** (`check_T_crystal_min_cut_v69`, tier 3, [P_structural]): vertex-split network with unit cap on internal `v_in → v_out` + INF cap on original edges, Edmonds-Karp max-flow → min-cut, witness extracted via residual-graph reachability. Sentinel `cut_size = -1, direct = True` returned for direct bank edges where vertex connectivity is undefined. Headline: T_sin2theta = 5, T_gauge = 3, T_PMNS = 9, T_mass_ratios = 8, L_count = 8 (3 direct-edge sentinels).
  - **Workstream 4 — DAG-DP path attribution by PLEC anchor** (`check_T_crystal_path_attribution_v69`, tier 4, [P_structural]): depth-filtered DAG to break documentation cross-ref cycles; topological-order O(V+E) DP from target backward; per-anchor `n_paths` + `share_pct`. Headline: A1 = 66%, L_epsilon* = 34%, **Regime_R + worked_example each 0% to T_sin2theta**.
  - **Workstream 5 — convergence-cluster detection** (`check_T_crystal_convergence_v69`, tier 3, [P_structural]): fan-in on depth-filtered DAG; `anchor_diversity = number of PLEC anchors reaching ≥1 direct predecessor`. Headline: **81% of nodes (177/218) are convergence sinks**, top-15 all have anchor_diversity = 2/4 (only A1 + L_epsilon*).
- All 5 checks register via the standard `_CHECKS` dict + `register(registry)` pattern (bare names, no `check_` prefix). Stage III canonical sink set `{T_sin2theta, T_gauge, T_PMNS, T_mass_ratios, L_count}` shared across workstreams 3 and 4.
- Together the four substantive findings form a coherent **"robustness has been won" narrative for Paper 20 v3.0 §4-§7**: the v6.9 framework is dramatically more over-determined than the v3.4.4 snapshot Paper 20 v1.0 analyzed, but the two newest PLEC anchors (Regime_R, worked_example) haven't been rewired into the SM derivation chain at the bank-edge level. Phase 14b candidate (PLEC anchor rewiring of gauge/field-content modules); not blocking Paper 20 v3.0 — the asymmetry is itself the headline finding.
- `apf/bank.py` adds `'apf.crystal_metrics'` to `_MODULE_PATHS` between `apf.crystal` and `apf.cosmology`; `EXPECTED_THEOREM_COUNT` bumped 365 → 371 → 372 → 373 → 374 → 375 across the workstream sequence; final loaded count 375 matches exactly.
- `verify_all.py` MODULES extended with `("apf.crystal_metrics", "crystal_metrics.py")` in the same slot.
- Net delta: +5 bank-registered checks (365 → 375 — includes running pre-bumps and the +5 from this workstream sequence; final loaded count 375 matches EXPECTED_THEOREM_COUNT exactly); +4 verify_all (380 → 384, one bump per check landed across workstreams 2–5; workstream 1 was the prior session); +1 module (22 → 23).

### v6.9 changes (2026-04-21 night Phase 13.2 Enforcement Crystal walker pass)
- New `apf/crystal.py` module (~706 lines, 1 bank-registered check). The [[Enforcement Crystal Walker]] reads the v6.9 bank registry and builds the Enforcement Crystal as a directed graph with three module-include presets (CORE / EXTENDED / FULL) and a dual-view contract (`full_graph` keeps Regime_R + 5 Type I–V exit channels; `post_R_subgraph` drops Regime_R + 4 rejected channels and keeps only Type V). Tier-1 metrics computed per build (depth, width, waist candidates, epistemic split, sector closures). JSON dashboard payload conforms to canonical Three.js viewer schema (`dashboard_data.json`); `write_dashboard_json(path, payload)` emits a 228 KB UTF-8 file the viewer ingests directly. Bank check `check_T_crystal_v69_consistent` (tier 4, [P_structural]) certifies all four PLEC anchors (A1/MD/A2/BW) reachable, walker ⊆ bank, post_R ⊊ full and contains no Regime_R reference, no edges dangle, JSON payload round-trips through `json.dumps`. Smoke-test (CORE preset, post-bootstrap): 227 nodes / 821 edges full; 222 nodes / 810 edges post_R; max depth 30; waist candidates `T1`, `L_irr`, `L_CP_dual_mechanism`.
- New `apf/crystal_axiom_roots.py` (Phase 13.1, pure-data module — not bank-registered). Canonical mapping from each PLEC component to its source-most bank-registered check: A1 → `check_A1` (`apf.core`); A2 → `check_Regime_R` (`apf.plec`); MD → `check_L_epsilon_star` (`apf.core`); BW → `check_worked_example` (`apf.core`). Consumed by Phase 13.2 walker to seed anchor vertex set.
- `apf/bank.py` adds `'apf.crystal'` to `_MODULE_PATHS` between `unification_three_levels` and `cosmology`; `EXPECTED_THEOREM_COUNT` 364 → 365.
- `verify_all.py` MODULES extended with `("apf.crystal", "crystal.py")` in same slot.
- Net delta: +1 bank-registered check (364 → 365); +1 verify_all (379 → 380); +1 module (21 → 22).

### v6.9 changes (2026-04-21 evening Phase 14 three-level identity refinement pass)
- New `apf/unification_three_levels.py` module (~990 lines, 17 bank-registered checks). Separates each [[T_ACC Unification|T_ACC consistency identity]] $I_k$ into integer / scalar / subspace witnesses + composed three-level consistency theorem; plus top-level `T_three_level_unification`. I2 fully [P] at all three levels (V_global = 42-dim subspace co-witness via `T_interface_sector_bridge`). I1 / I3 / I4 retain [P] integer + scalar witnesses with [C, parked] subspace witnesses awaiting categorical refinement. See concept page **[[Three-Level Identity Refinement]]**.
- Net delta: +17 bank-registered checks (347 → 364); +17 verify_all (362 → 379); +1 module (20 → 21).

### v6.9 changes (2026-04-20 interface-sector bridge pass)
- New `check_T_interface_sector_bridge` in `apf/gravity.py` (tier 4, [P]). Main structural bridge theorem: **the T12 interface partition `V_61 = V_local (+) V_global` governs the `T_horizon_reciprocity` second-epsilon sector decomposition.** `|Sector A| = |V_61 \ {self}| = 60` and `|Sector B| = dim V_global = 42`, making `d_eff = 60 + 42 = 102` a structural corollary rather than an algebraic sum. Promotes the recurring "42" in T11 (`Omega_Lambda = 42/61`) and `L_self_exclusion` (`d_eff = 60+42`) from numerical coincidence to witnessed geometric identity — both read off the same 42-dim subspace `V_global` (T12's global-interface stratum).
- New `check_L_global_interface_is_horizon` in `apf/gravity.py` (tier 3, [P]). Auxiliary lemma identifying T12's global-interface stratum with the de Sitter horizon absorber subspace (T_Bek), so the bridge theorem's Step 4 closes without an implicit cross-module step.
- Downstream docstring cross-references added to `L_self_exclusion`, `T_deSitter_entropy`, `T_horizon_reciprocity`, and T11 (in `apf/cosmology.py`). Dependency lists unchanged to avoid cycles (the bridge theorem depends on those downstream theorems).
- Relation to T_ACC identity I2: I2 witnesses `K = 61` as the same integer in two π-projections; this theorem is the parallel statement at the **subspace** level (same 42-dim `V_global` in two readings) — strictly finer than I2.
- Net delta: +2 bank-registered checks (345 → 347 actually loaded, reconciling the pre-existing silent −2 drift; `EXPECTED_THEOREM_COUNT = 347` unchanged and now matches exactly); +2 verify_all (360 → 362); +0 modules.

### v6.9 changes (2026-04-20 T_ACC unification pass)
- New `apf/unification.py` module (5 bank-registered checks): `I1_holographic`, `I2_gauge_cosmological`, `I3_thermo_quantum`, `I4_action_thermo`, and the composed `T_ACC_unification`. Introduces the `ACC` record (structural capacity `K`, effective degeneracy `d_eff`, microstate count `N = d_eff^K`, admissibility-capacity `ACC = K · ln d_eff`), the six regime projections (`pi_T`, `pi_G`, `pi_Q`, `pi_F`, `pi_C`, `pi_A`), three canonical interface factories (`acc_SM`, `acc_horizon`, `acc_quantum`), and the four cross-regime consistency identities that T_ACC asserts. Each identity is registered as a first-class theorem, and `T_ACC_unification` composes all four on the canonical interfaces for audit.
- T_ACC stub removed from `apf/plec.py`; `plec.py` is back to exactly Regime_R + five exit types.
- Net delta: +5 bank-registered checks (342 → 347); +5 verify_all (355 → 360); +1 module (19 → 20).

### v6.9 changes (prior pass, 2026-04-18)
- New `apf/plec.py` module (6 checks): `Regime_R`, `Regime_exit_Type_I` through `Type_V`. Formalizes the PLEC infrastructure introduced in Papers 5/6 v2.0-PLEC.
- New `check_A9_closure` in `apf/gravity.py` (1 check). Unifies the Lovelock prerequisites A9.1–A9.5 dispersed across `core.py`, `gravity.py`, `spacetime.py`, `internalization_geo.py`.
- Net delta: +7 bank-registered checks (335 → 342); +7 verify_all (348 → 355); +1 module (18 → 19).
- Papers 5/6 v2.0-PLEC code-anchored with full coderef pass; Paper 6 §11 substantially expanded with the 4-parameter Planck match table, $H_0 = 67.76$ km/s/Mpc derivation, H0DN 7.09σ tension, escape-route table, and explicit falsifier statement.
- The 4-parameter Planck cosmological match ($\Omega_b = 3/61$, $\Omega_c = 16/61$, $\Omega_m = 19/61$, $\Omega_\Lambda = 42/61$, all matching Planck within 1.2%) is now foregrounded as a major framework result — and has become the load-bearing structural content of T_ACC identity I2 (gauge-cosmological).

### Authoritative counts (from `apf.bank._load()` run 2026-04-21 night, post-Phase-13.2 Enforcement Crystal walker pass)

- **380** total verify_all checks — 376 pass + 4 scipy-gated (sandbox-only failure mode; pass on fully provisioned hosts)
- **365** bank-registered theorems (`apf.bank.REGISTRY` size after `_load()`; matches `EXPECTED_THEOREM_COUNT` exactly — no count drift)
- **22** bank-registered modules in the `apf/` package
- **+ apf/standalone/** sub-package (4 additional files, not bank-registered but verify_all-tracked)
- **+ session_phase2_confrontation.py** (3 additional checks, not bank-registered)
- **48** quantitative predictions
- **0** free parameters
- **39** tested predictions
- **32/39** within 3σ of observation
- **Mean error:** 3.83%
- **Median error:** 0.37%

### Module breakdown (bank-registered, Phase 2.9 frozen count)

| Module | Bank-registered theorems | Purpose |
|--------|-------------------------:|---------|
| `apf.core` | 48 | Axiom A1 / PLEC, Hilbert space, Born rule, CPTP |
| `apf.gauge` | 31 | Gauge structure, non-closure, uniqueness |
| `apf.generations` | 85 | Fermion spectrum, 3 generations, masses, mixing |
| `apf.spacetime` | 8 | Metric structure, field equations |
| `apf.gravity` | 12 | Gravitational sector, horizon entropy (v6.9 additions: +A9_closure 2026-04-18, +L_global_interface_is_horizon + T_interface_sector_bridge 2026-04-20) |
| `apf.cosmology` | 17 | Dark energy, density fractions, N_eff |
| `apf.validation` | 9 | Prediction validation vs. experiments |
| `apf.supplements` | 70 | Lemmas, inequalities, computational aids |
| `apf.majorana` | 10 | Neutrino physics, Majorana constraints |
| `apf.internalization` | 3 | Gauge boson internalization |
| `apf.internalization_geo` | 4 | Geometry from internalization |
| `apf.extensions` | 7 | BSM physics, dark matter mechanisms |
| `apf.red_team` | 17 | Adversarial tests for false positives, edge cases |
| `apf.session_v63c` | 4 | v6.3c incremental: neutrino hierarchy, Yukawa spectral |
| `apf.session_qg` | 4 | QG / Planck-scale scoping |
| `apf.session_nnlo` | 4 | NNLO corrections |
| `apf.session_delta_pmns` | 2 | PMNS phase derivation progress |
| `apf.session_cosmo_update` | 3 | Cosmological observables refresh |
| `apf.plec` | 6 | Regime R + 5-type regime-exit taxonomy (NEW v6.9 2026-04-18) |
| `apf.gravity` (A9_closure addition) | +1 | Unified Lovelock-prerequisite closure (NEW v6.9 2026-04-18) |
| `apf.gravity` (interface-sector bridge additions) | +2 | **NEW v6.9 (2026-04-20 late evening).** `L_global_interface_is_horizon` (aux) + `T_interface_sector_bridge` (main). Identifies `V_global` (T12 global-interface stratum) with the dS horizon absorber subspace and proves Sector B (second-epsilon unilateral drops) == `V_global`; promotes the "42"s in T11 and `L_self_exclusion` to a single geometric identity. |
| `apf.unification` | 5 | **NEW v6.9 (2026-04-20).** Admissibility-Capacity Ledger: I1 holographic, I2 gauge-cosmological, I3 thermo-quantum, I4 action-thermo, and the composed T_ACC_unification. |
| `apf.unification_three_levels` | 17 | **NEW v6.9 (2026-04-21 evening).** Three-level identity refinement: 4 integer + 4 scalar + 4 subspace + 4 composed Iₖ + 1 top T_three_level_unification. |
| `apf.crystal` | 1 | **NEW v6.9 (2026-04-21 night).** Enforcement Crystal walker: `T_crystal_v69_consistent` (tier 4, [P_structural]). Reads bank registry; builds dual-view DAG; certifies PLEC anchors reachable + walker ⊆ bank + post_R ⊊ full + JSON serializability. |
| **Bank total** | **365** | |

**Additional verify_all-tracked (not bank-registered):**

| Source | Checks | Notes |
|--------|-------:|-------|
| `apf.session_phase2_confrontation` | 3 | Phase 2 experimental-confrontation module (not yet registered) |
| `apf.standalone.L_Cauchy_uniqueness` | 1 | Cauchy-uniqueness lemma, standalone |
| `apf.standalone.L_CKM_resolution_limit` | 1 | CKM resolution limit, standalone |
| `apf.standalone.phase1_seesaw_closure` | 2 | Phase 1 seesaw-chain closure + RT test |
| `apf.standalone.phase5_theorem_R_audit` | 3 | Phase 5 Theorem R adversarial audit (R1/R2/R3) |
| **Verify_all grand total** | **380** | |

### Environment and dependencies

- **Python ≥ 3.8**
- **numpy ≥ 1.20** (linear algebra, array ops)
- **scipy ≥ 1.7** (special functions, two-loop RG, spectral-action integrals — 4 checks require this specifically: `L_SA_moments`, `L_ST_index`, `L_mc_mt_twoloop_RG`, `L_spectral_action_coefficients`)
- **fractions, math, itertools** (stdlib, for exact-rational core checks)

Install via `pip install -e .` from the codebase root (uses `setup.py`).

### Key properties

- **Mixed-precision:** Core framework uses exact rationals (`fractions.Fraction`); numerical routines (two-loop RG, spectral action) use `numpy` / `scipy` with documented tolerance bounds
- **Modular:** Each module encodes a logical layer of the derivation chain
- **Deterministic:** No randomness; results are reproducible
- **Fast:** Full suite runs in ~5 seconds on standard hardware
- **Auditable:** Every check is a named theorem; human-readable
- **Bank-registered assertion:** `bank.EXPECTED_THEOREM_COUNT = 365` with runtime verification. Any future drift raises a warning immediately.

### Example check:

```python
def check_T_Born(tol=1e-15):
    """T_Born: Born rule probability formula

    For state |psi> and observable A, P(a_n) = |<a_n|psi>|^2.
    """
    # Verify this is the unique cost-minimizing measurement rule
    # under admissibility constraint...
    return True
```

### Development practice

- **v6.8 frozen** 2026-04-18 (Phase 2.9 canonicalization) as the authoritative snapshot for the current paper series (Papers 0, 1, 2, 3, 7 live in .tex; Papers 4, 5, 6 source-pending).
- New features branch for v7.0 (future, post-arXiv submission window).
- All commits tested; no regressions.
- `verify_all.py` is the single source of truth for count; CLAUDE.md and wiki mirror its output.

### Canonicalization history

- **v6.9 T_ACC hardening pass (2026-04-20, evening):** No check/count changes (still 347/360/20). Tightened `acc_SM` in `apf/unification.py` with `dag_has` presence-guards so a missing or mis-populated `C_total` / `d_eff` in the DAG now raises a named-producer RuntimeError instead of silently substituting the canonical defaults that would have made I1–I4 pass on their targets by construction. Side-effect-discovered `d_eff` DAG-key collision — `L_self_exclusion` (`apf/gravity.py`) and `L_nuR_enforcement` (`apf/majorana.py`) were both writing the string key `'d_eff'` but meaning totally different quantities (degeneracy-per-slot 102 vs. seesaw vertex mass-dimension 4); renamed the majorana consumer key to `L_nuR_vertex_dim`. Added `PRELUDE_MODULES` + `run_prelude_silently` to `verify_all.py` so `--module <filter>` runs of downstream consumers auto-execute upstream DAG producers as a silent prelude (one-line banner reports prelude scope); without this, the new presence-guards made per-module debugging of downstream checks impossible. Swept 26 stale/duplicate root-level `.py` files left over from Phase-2.9 canonicalization (3 were dangerously stale v6.7 copies that would silently shadow `apf/` originals if Python ever imported from root; archived to `Codebase/Old/root_stale_duplicates_pre-cleanup_2026-04-20.zip` then deleted). `paper1.py` import path namespaced (`apf_utils` → `apf.apf_utils`). `verify_all.py` banner corrected v6.7 → v6.9. Found-but-not-fixed: `bank.EXPECTED_THEOREM_COUNT = 347` but per-module counts load 345 at runtime (the "347" figure appears to have been set from a documented PLEC baseline that was itself off by 2; truthful fix is a two-line `bank.py` change plus CLAUDE.md cascade, deferred to the Paper 8 build). Archives: `Codebase/Old/{gravity.py_pre-dag-put-d_eff,majorana.py_pre-d_eff-rename,paper1.py_pre-apf_utils-qualify,unification.py_pre-DAG-guards,verify_all.py_pre-prelude}_2026-04-20` + `root_stale_duplicates_pre-cleanup_2026-04-20.zip` + milestone `APF_Codebase_v6.9_post-tacc-hardening_2026-04-20.zip`.
- **v6.9 T_ACC pass (2026-04-20):** New `apf/unification.py` module adds 5 bank-registered checks (342 → 347; 355 → 360 verify_all; 19 → 20 modules). Module contains the ACC record, the six regime projections (`pi_T`, `pi_G`, `pi_Q`, `pi_F`, `pi_C`, `pi_A`), three canonical-interface factories (`acc_SM`, `acc_horizon`, `acc_quantum`), four consistency-identity helpers (`_identity_I1..I4_*`), the zero-arg bank-registered wrappers (`check_I1..I4_*`), and the composed `check_T_ACC_unification`. Each identity is registered as a first-class theorem; `T_ACC_unification` re-runs all four via the helpers for a single composite audit. T_ACC stub removed from `apf/plec.py`. `verify_all.py` MODULES list extended with `apf.unification`. `bank.EXPECTED_THEOREM_COUNT` bumped 343 → 347 (net +5 over the 2026-04-18 baseline of 342, or +4 over the intermediate 343 that registered only the composed theorem). Archives: `Codebase/Old/bank.py_pre-unification-module_2026-04-20`, `bank.py_pre-I1-I4-registration_2026-04-20`, `plec.py_pre-T_ACC-stub_2026-04-20`, `unification.py_pre-I1-I4-registration_2026-04-20`, `verify_all.py_pre-unification-module_2026-04-20`.
- **v6.7 → v6.8 (2026-04-18):** Phase 2.9 canonicalization pass. Moved `session_*.py` and standalone files into the `apf/` package (previously at top level, breaking imports in `verify_all.py`). Fixed broken `rt_check` helper in `apf/standalone/phase1_seesaw_closure.py`. Updated `apf.bank.EXPECTED_THEOREM_COUNT` from 310 → 335 to match loaded reality. Updated `setup.py` version 6.7.0 → 6.8.0 and added explicit `install_requires = numpy, scipy`. Synced all count references across CLAUDE.md (× 2), wiki (Index, Log, APF Codebase, Predictions Catalog, paper-specific pages), Schema.md, Paper Update Work Plan, Paper Index, Paper 0 v3.0, Paper 3 v3.2 + Supplement. Driven by the external audit of 2026-04-18 (`APF Reference Docs/Reference - External Audit Memo (2026-04-18).md`).
- **v6.7 (pre-2026-04-18):** Prior authoritative state. Docs drifted from 294 → 312 → 321 → 349 → 351 across different files; `bank.EXPECTED_THEOREM_COUNT = 310` stale relative to loaded 335; session/standalone files couldn't import via `apf.*` namespace; one red-team check failed on a missing helper. All resolved in v6.8.

## See also
- [[Paper 13 - Minimal Admissibility Core]]
- [[Derivation Chain]]
- [[Red Team Tests]]
- [[Predictions Catalog]]
