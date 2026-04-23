---
type: paper
domain: apf
layer: cross-cutting
created: 2026-04-14
updated: 2026-04-22 (Phase 13–14c master-reference update: v8.4 → v8.5)
sources:
  - Paper_13_Minimal_Admissibility_Core_v8.5.tex
  - Edit Proposals - Paper 13 Descriptive Framing.md
---

# Paper 13 - Minimal Admissibility Core
## "The Master Reference"

**Status:** **v8.5 (69pp, Phase 13–14c master-reference update, 2026-04-22).** Codebase anchors refreshed to v6.9 post-Phase-14c (376 bank-registered theorems across 23 modules + standalone sub-package, 385 total verify_all checks, 48 quantitative predictions, 0 free parameters) throughout abstract / front-matter scorecard / derivation-chain prose / DAG prose / architecture discussion / migration note / crystal-metrics intro. Module Architecture table gains four new rows for the Enforcement Crystal mechanization: `crystal_axiom_roots.py` (Phase 13.1, PLEC anchor constants), `crystal.py` (Phase 13.2, bank walker with dual-view contract), `crystal_metrics.py` (Phases 13.3 + 14c, five bank-registered analytical-metrics checks including SCC-condensed path attribution), and `unification_three_levels.py` (Phase 14, 17 integer/scalar/subspace witnesses of I1–I4). Theorem Reference Index gains three matching sub-tables. Appendix J version history extended with v8.5 row documenting Phase 13–14c delta + Phase 14b rewiring + Phase 14c SCC backend. Previous: v8.4 (Phase 12.6, 2026-04-21).

## Changelog

- **v8.5 (2026-04-22):** Phase 13–14c master-reference update. Four structural additions: (1) front-matter scorecard 347/362/20 → 376/385/23 (48 predictions unchanged); date-line extended with "Phase 13–14c Enforcement Crystal mechanization: April 2026." (2) Module Architecture table gains four rows: `crystal_axiom_roots.py` (PLEC anchor constants A1/MD/A2/BW → source-most checks, data-only, Phase 13.1), `crystal.py` (bank-to-graph walker + dual-view contract + three module presets, Phase 13.2), `crystal_metrics.py` (Brandes betweenness + single-deletion cascade + Menger min vertex cuts + SCC-condensed DAG path attribution via Tarjan 1972 + convergence-cluster detection, Phases 13.3 + 14c), and `unification_three_levels.py` (17 integer/scalar/subspace witnesses of I1–I4 + composed three-level identities + $T_{\mathrm{three\,level\,unification}}$, Phase 14). (3) Theorem Reference Index adds three new sub-tables matching the new modules; `unification_three_levels.py` sub-table notes $I2$ fully [P] at all three levels via $V_{\mathrm{global}}$ 42-dim co-witness; `crystal_metrics.py` sub-table notes the 22-node SM-core SCC with pointer to Paper 20 v3.1 Finding 5. (4) Version History entry for v8.5 documents Phase 14b (Regime_R + worked_example as bank-edge predecessors on SM spine — collapses Paper 20 v3.0 "zero paths to $\sin^2\theta_W$" asymmetry by construction) and Phase 14c (Tarjan SCC condensation as primary cycle-aware path-attribution backend). xelatex × 3 clean. 69pp (from 67pp). Archive: `Old/Paper_13_Minimal_Admissibility_Core_v8.4_pre-v8.5.{tex,pdf}`. Pre-existing `§??` from unresolved `\ref{sec:field_content}` inherited from v8.4 — not a v8.5 regression, flagged for a later correctness pass.
- **v8.4 (2026-04-21):** Phase 12.6 master-reference update. Five targeted edits: (1) front-matter scorecard 129/8/827 → 347/362/20-mod/48 predictions; (2) Module Architecture table expanded to 20 modules + standalone sub-package with entries for the two v6.9 additions `plec.py` (Regime R + five-type regime-exit taxonomy, 2026-04-18) and `unification.py` (ACC ledger, 2026-04-20); (3) Theorem Reference Index gains new sub-tables for `unification.py` (I1 holographic, I2 gauge-cosmological, I3 thermo-quantum, I4 action-thermo, T_ACC_unification) and `plec.py` (Regime R, Types I–V); `gravity.py` sub-table extended with 6 new rows (bridge theorem + auxiliary lemma + T_Bek + T_horizon_reciprocity + T_graviton + A9_closure); (4) scattered stale prose counts refreshed throughout body; (5) Version History entry for v8.4. xelatex × 3 clean. 67pp (from 64pp). Archive: `Old/Paper_13_Minimal_Admissibility_Core_v8.3-PLEC_pre-Phase12.6.{tex,pdf}`. One pre-existing undefined `\ref{sec:field_content}` warning inherited from v8.3-PLEC (not introduced).
- **v8.3-PLEC (2026-04-17):** Phase 2.8a surgical PLEC propagation. Header/fancyhead/title retagged v8.3-PLEC; epigraph (already the canonical PLEC lay statement) annotated with a framing note identifying PLEC and its four structurally necessary components (A1, MD, A2, BW); abstract expanded with a PLEC-alignment sentence pointing to where the non-A1 components enter (L_col argmin, C_total=61 count); L_col subsection gains a "note on PLEC alignment" reconciling Paper 13's "derived from A1" shorthand with Paper 1 v4.0-PLEC's treatment of A2 as structurally necessary; v8.3-PLEC row added to Appendix J version history. No theorems, numerics, or code references changed. Full v9.0 restructure around PLEC's four-component inventory deferred. Compiled clean via xelatex (three passes), 64pp.
- **v8.2 (2026-04-16):** Descriptive-framing pass applied. 18 edits (P13.0–P13.17) from [[Edit Proposals - Paper 13 Descriptive Framing]]: new abstract reader's note, new `\subsection{Descriptive versus selective readings}` framing box (§2.4), section title rewrites ("What the Axiom Forces" → "What the Axiom Admits"; "Minimality — nature is parsimonious" → "Minimality — the admissibility set is an $\arg\min$"), and fourteen running-prose rewrites converting "selects/forces/selection" to "admits/unique admissible/admissibility". Incidental fix to pre-existing `\not\xrightarrow` syntax in §A.3 (rewritten as prose). Zero formal content changes. Compiled via xelatex (paper uses fontspec); two passes clean; 63pp.
- **v8.1 (pre-April 2026):** Comprehensive self-contained exposition. 18 modules, 335 bank-registered theorems, 348 verify_all checks (counts canonicalized to codebase v6.8 on 2026-04-18; v8.1 was originally quoted as "294+ theorems, 349 checks" which turned out to be stale metadata against the actual loaded count).

**Core claim:** A comprehensive, self-contained exposition of the entire APF framework in minimal form. All 23 bank-registered modules + apf/standalone/, 376 bank-registered theorems, 385 verify_all checks (codebase v6.9: 2026-04-22 Phase 14c SCC-condensed path attribution + Phase 14b PLEC anchor rewiring + 2026-04-21 Phase 13.1–13.3 Enforcement Crystal mechanization + Phase 14 three-level identity refinement + 2026-04-20 interface-sector bridge + 2026-04-20 T_ACC unification + 2026-04-18 PLEC formalization), unified derivation from [[Axiom A1]]. The canonical reference for understanding the full architecture.

**Structure:**
- Layer 1 (Spine): Axiom, operability, canonical object
- Layer 2 (Structure): Non-closure, gauge uniqueness
- Layer 3 (Ledgers): Cost ledgers, irreversibility (abbreviated)
- Layer 4 (Constraints): Field content, fermion spectrum, masses
- Layer 5 (Quantum): Hilbert space, Born rule, CPTP evolution
- Layer 6 (Dynamics): Spacetime, gravity, cosmology
- Layer 7 (Action): Internalization, Lagrangian, BSM
- Cross-cutting: Predictions, open problems, reconstruction pathways

**Content:**
- Unified notation and conventions
- All core theorems with proof sketches
- Module-by-module codebase map
- 47+ predictions with experimental status
- Open problems and future directions
- Reading guide for papers 1-7

**Implementation:** All 18 modules reference v8.1

**Status notes:**
- Comprehensive and stable
- Acts as a master index
- Suitable as a standalone reference

## See also
- [[Derivation Chain]]
- [[Papers 1-7|Paper 1 - Spine]]
- [[APF Codebase]]
- [[T_ACC Unification]] — `unification.py` + `unification_three_levels.py` sub-tables in v8.5 Theorem Reference Index
- [[Interface-Sector Bridge]] — `gravity.py` additions in v8.4 Theorem Reference Index (carried forward in v8.5)
- [[Paper 20 - The Enforcement Crystal]] — v3.1 is the companion paper for the `crystal.py` + `crystal_metrics.py` + `crystal_axiom_roots.py` mechanization documented in v8.5 Appendix E
