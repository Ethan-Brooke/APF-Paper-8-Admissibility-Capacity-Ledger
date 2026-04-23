---
type: paper
domain: apf
layer: 0-ontology
created: 2026-04-14
updated: 2026-04-21 (v3.0 → v3.1 — three-layer-reading pass: new §8 subsection "The Three Readings of an Interface Inside Regime R" [architecture K=61 / commitment 2ε per slot / target d_eff=102 partitioned 60+42]; orthogonality paragraph in §3.4 vs. four ontological layers; new §5 subsection "What Regime R makes possible: rigidity of architecture and commitment"; new §sec:T_ACC_commitments item registering V_global = Sector B = horizon absorber as strictly finer than I2; Book II projection additions for Papers 2/5/6/7 naming which layer each paper's projection reads)
sources: []
---

# What Physics Permits
## Paper 0 — An Ontology for Admissibility Physics

**Status:** **v3.1 (2026-04-21), 103pp .tex/.pdf** (up from 98pp at v3.0 close of 2026-04-20). v3.1 is the first paper-side landing of Phase 12's three-layer-reading propagation; content and proof structure from v3.0 unchanged, six structural additions aligned to the codebase v6.9 bridge theorem + T_ACC unification + upstream reference doc `Reference - Three-Layer Ontology and 61-102 Dynamical Picture (2026-04-21).md`. v3.0 was the full-structure rebuild (2026-04-17) replacing the Feb-2026 v2.0 .docx "gateway" narrative; substantially expanded over the April 19–20 sessions with three new chapters, three new figures, Book II chapter fills, T_ACC ledger chapter, falsifier aggregation upgrade, and ~30 new bibliography entries. v3.1 is the ontology paper of the APF series: the reference statement of what the framework commits to, from which all technical papers derive.

**v3.1 changes (2026-04-21):**
- New §8 (Admissibility-Capacity Ledger) subsection `sec:architecture_commitment_target` — "The Three Readings of an Interface Inside Regime R," a load-bearing apfbox (rewritten to plain-English prose 2026-04-21 after an initial question-based draft was rejected as stylistically off) defining the architectural layer (K = 61 channel-slot count, regime-rigid under L_count), the commitment layer (2ε per slot per instant, regime-forced under L_nc + L_irr), and the target layer (d_eff = 102 per slot, partitioned into Sector A [60, bilateral anchors] + Sector B [42, V_global / horizon absorber] by the bridge theorem).
- §3.4 (four-layer hierarchy) orthogonality paragraph distinguishing the three-reading insert from the existing four ontological layers (keystone / A1 / MD-BW / A2).
- §5 (Regime R) new subsection `sec:regime_R_rigidity` — "What Regime R makes possible: rigidity of architecture and commitment" — frames R's positive content in three-layer terms, with each of the five exit types mapped to the rigidity condition it dissolves.
- §sec:T_ACC_commitments — new item 3 registering the subspace-level identity V_global = Sector B = horizon absorber as strictly finer than T_ACC identity I2 (integer-level K = 61). Cites `check_T_interface_sector_bridge` + `check_L_global_interface_is_horizon`.
- Book II Papers 2 / 5 / 6 / 7 projection sections each extended with a paragraph naming which layer the paper's projection reads: Paper 2 primary-architectural + target-partition supplier, Paper 5 target layer (quantum apparatus; measurement as Sector A↔Sector B crossing), Paper 6 joint architectural+target (Ω_Λ = 42/61 as subspace fingerprint), Paper 7 commitment layer (2ε per slot per instant, 122 ε·dt universe-instantaneous identity, Z_0 factorised into all three readings).
- T_ACC table inside §sec:T_ACC_statement switched from plain tabular to tabularx with footnotesize; fixes pre-existing overflow from the apfbox frame.
- Version bump v3.0 → v3.1 in preamble comment + title block; detailed changelog in preamble.
- Pdflatex × 3 passes clean at 103pp.

**Core claim:** The framework's ontological commitments are a four-layer hierarchy: (1) the keystone — *meaning requires enforceable distinction*; (2) finiteness of enforcement (A1); (3) admissibility structure (A1 + MD + BW); (4) realization as the $\arg\min$ over the admissible set (A2; the [[Principle of Least Enforcement Cost]] unifies all four). PLEC is stated canonically: *reality is the minimum-cost expression of distinction compatible with finite capacity*. This sentence is the framework's frozen lay handle.

**Three-book structure (updated 2026-04-19):**

**Book I — The Ontology** (Ch 1–7, expanded from 5). Substantive core.
- Ch 1 The Ontological Keystone. Meaning requires enforceable distinction. Why the keystone is prior to dynamics and prior to representation. What the keystone rules out (zero-cost distinctions, distinctions with no consequence, notational fictions). **Now includes explicit structural-realist anchor (Worrall, Ladyman–Ross, French).**
- Ch 2 The Finiteness of Enforcement (A1). A1 stated with its capacity-bound formula. What A1 does (admissible set from above). What A1 does not do. Stanford-reviewer countermodel separating A1 from MD. **A1-support paragraph softened from "establish" to "motivate a transfer" framing; new citations to Bousso, Wheeler, Adlam, Spencer-Brown.**
- Ch 3 The Admissibility Structure. A(ρ, Γ) with all four bounds. MD, BW. Four-layer hierarchy. Structural independence of all four components.
- Ch 4 The Principle of Least Enforcement Cost (PLEC). Formal statement with canonical lay sentence. The argmin is a locator, not a process. Scope of the argmin. PLEC in the variational tradition. Non-merge discipline. **Now anchored to Curiel, Ruyant–Guay, Butterfield on Lagrangian possibilities + Lange/Adlam on laws-as-constraints.**
- **Ch 5 Regime R and the Regime-Exit Taxonomy (NEW 2026-04-19).** Regime R's four conditions R1–R4. The five-type exit taxonomy (Type I collapse; Type II minimizer nonuniqueness; Type III cross-class; Type IV loss of smoothness/locality; Type V pure representational redundancy — with three-flavor distinction per Healey / Nguyen–Teh–Wells). Exhaustiveness-and-no-overlap. How the taxonomy plays across the series. Includes figure.
- Ch 6 Descriptive Reading (was Ch 5). The convention stated. Applied glossary of narrative↔descriptive translations.
- **Ch 7 Philosophical Positioning (NEW 2026-04-19).** Three readings of PLEC (meta-law / reconstruction / ontology) mapped onto the four-layer hierarchy; explicit ordering commitment ("ontology first, reconstruction generative, meta-law explanatory"). Bridge principles for cross-domain transferability. Engagement with principal objections (Brown–Timpson, Hagar–Hemmo, Timpson, Felline, Healey, Nguyen–Teh–Wells, Goyal, Berghofer).

**Book II — Derivational Commitments** (Ch 6–13). One chapter per paper under fixed six-section template (Scope / Commitments drawn on / Structural results / Where "forces" is a descriptor / What the paper does not claim / Status). (Chapter numbering has shifted: old 6/7/8 now 8/9/10 etc. due to Ch 5 insertion; effective Book II slots for Papers 1, 2, 3, 4, 5, 6, 7, 8 follow the convention.)
- Paper 1 (Enforceability of Distinction) — full.
- Paper 2 (Structure of Admissible Physics) — full.
- Paper 3 (Ledgers) — full.
- **Paper 4 (Constraints — Field Content) — EXPANDED 2026-04-19** from placeholder. Full six-section template at v2.0 resolution: 19 structural results covering gauge template uniqueness, 1-of-4680, 45-fermion spectrum, sin²θ_W = 3/13, Cabibbo, CKM, PMNS, dark matter, cosmological constant, 5 falsifiers. Four paper-specific hypotheses H1–H4.
- **Paper 5 (Quantum Structure) — EXPANDED 2026-04-19** from placeholder. 9 structural results covering linearity, complex Hilbert via Frobenius trichotomy, T_tensor, T_Born via Gleason, T_Hermitian, T_CPTP with Stinespring interface-mediator reading, T_decoherence as Type III exit, L_spectral_action_internal, L_QG_P1_closure UV completion. H1/H2/H3 hypotheses.
- **Paper 6 (Dynamics and Geometry) — EXPANDED 2026-04-19** from placeholder. 20 structural results covering PLEC path form, regime-exit taxonomy, correlation-cost geodesics, T_7B, T_8 via six-Δ chain, A9 closure + Lovelock + Einstein equations, T_11 Ω_Λ = 42/61, four-parameter Planck match (with figure), H0 tension falsifier, T_12/T_12E dark matter, BH thermodynamics. H1–H4 hypotheses.
- Paper 7 (Action, Internalization, and the Lagrangian) — full.
- **Paper 8 (The Admissibility-Capacity Ledger) — EXPANDED 2026-04-20** from placeholder. Full six-section template. Scope: the framework's unification paper — every regime-specific cost or capacity quantity (entropy, horizon area, gauge capacity, cosmological fractions, action partition) is a projection of a single framework-level admissibility-capacity. Commitments: two-scalar structure ($K$ structural capacity, $\mathrm{ACC} = K \ln d_{\mathrm{eff}}$ log-count); six regime projections ($\pi_T$, $\pi_G$, $\pi_Q$, $\pi_F$, $\pi_C$, $\pi_A$); four consistency identities (holographic, gauge-cosmological, thermo-quantum, action-thermo). Structural results: T_ACC as the unification theorem; at the SM, $K_{SM} = 61 = 45 + 4 + 12$, $d_{\mathrm{eff}} = 102$, $\mathrm{ACC}_{SM} \approx 282.12$ nats = $S_{dS}$. Falsifiers: projection-map failure (any framework-derived quantity that refuses the ACC form); denominator failure (H0DN 7.09σ is the current active test). Registered at codebase level via `apf/unification.py` with five bank-registered checks (I1/I2/I3/I4 + composed T_ACC_unification); the full formal paper is reserved for the Paper 8 build.

**Book III — The Framework in the World** (Ch 14–17, expanded from 14–16).
- **Ch 14 Relation to Prior Work (NEW 2026-04-19).** Principle-theory tradition (Einstein, Brown–Timpson, Felline). Information-thermodynamics line (Jaynes, Landauer, Bekenstein, Bousso, Wheeler, Timpson). Quantum reconstruction — now split into technical reconstructions (Hardy, CBH, CDP, Dakić–Brukner) and philosophical reflection on reconstruction (Goyal, Berghofer, Felline). Variational principles and meta-laws (Lange, Adlam, Butterfield, Curiel, Ruyant–Guay). Gauge redundancy and structural ontology (Healey, Nguyen–Teh–Wells, Ladyman–Ross). Further reading by entry point.
- Ch 15 Where Standard Physics Lives (was Ch 14). Regime scoping.
- Ch 16 Status, Predictions, Falsifiers (was Ch 15). Codebase v6.9 scorecard (347 theorems, 360 checks, 48 predictions, 32/39 within 3σ, 3.83% mean, 0.37% median). Representative predictions table. **Falsifier aggregation upgraded 2026-04-19: H0 tension (7.09σ, Routes I–IV excluded, Route V insufficient by ~2×); PMNS-T2K F6; four-parameter Planck-match as 4-axis falsifier.** **2026-04-20: four-parameter Planck match promoted to the load-bearing content of T_ACC identity I2 (gauge-cosmological); the 7.09σ H0DN tension is now a denominator-level falsifier against T_ACC, not only Paper 6 §11.4.**
- Ch 17 How to Read the Papers (was Ch 16).

**Canonical sentence placements (load-bearing):**
- Title-page epigraph
- Book I Ch 4.2 (inside the PLEC canonical-statement tcolorbox)
- Book III Ch 16.6 (closing)

**What changed vs. v2.0:**
- v2.0 was a "gateway, not a textbook" five-act narrative (Preface + Act I The Problem + Act II The Reframe + Act III The Map + Act IV Where Standard Physics Lives + Act V The Reading Guide + Appendix Reference Table).
- v3.0 is the ontology paper, not a gateway. Four-layer hierarchy stated explicitly. PLEC named (v2.0 never used the word). Canonical sentence placed (v2.0 did not carry it). Paper 7 at v2.0 scope (spectral action / SM Lagrangian), not v1.0 scope (minimum action only). Paper 1 correctly described as carrying the full quantum reconstruction (v2.0 had said "does NOT derive any specific physics," now false under v4.0-PLEC). MD kept structurally separate from A1 (v2.0 elided this). Descriptive-framing discipline applied sentence-level, not via a global flag. Three-book structure replaces the five-act container so Book II can expand as technical papers mature.

**Audit driver (2026-04-17):** Ran full ontological audit against Paper 1 v4.0-PLEC + Paper 2 v5.3-PLEC + both Technical Supplements + Axiom Inventory v2.1. Found 9 blocker-level divergences (PLEC absence, keystone/A1 conflation, Paper 7 scope, Paper 1 scope, MD-as-separate, canonical sentence, five-vs-six layer inconsistency, etc.), 12 major, 4 minor. Ethan directed a full rebuild at v3.0 scope rather than a surgical port — "this is the ontology paper, words matter a lot, no shortcuts, take the time to craft precise language based on the math."

**Status notes:**
- v3.0 is shippable. Book II Chs 9–11 and 13 carry `pendingbox`-styled placeholders that will fill in as Papers 4, 5, 6, 8 mature; the structure is forward-compatible.
- pdflatex × 3 passes clean. No undefined references, no errors, 12 cosmetic overfull hboxes (description-list labels with long code refs, Paper 1 style inheritance).

**Archive:** `Old/Reference - What Physics Permits v2 Complete_pre-v3.docx` (the Feb-2026 .docx v2.0 that v3.0 replaces).

**Source:** `Papers/Paper 00 - What Physics Permits/Paper_0_What_Physics_Permits_v3.1.tex` + `.pdf`.

**Zenodo DOI (v2.0 deposit):** https://doi.org/10.5281/zenodo.18605692. A v3.0 deposit will follow once Ethan reviews the compiled PDF.

## See also
- [[Paper 1 - Spine]]
- [[Paper 13 - Minimal Admissibility Core]]
- [[PLEC Rollout Plan]]
- [[Derivation Chain]]
- [[Axiom A1]]
- [[Principle - Minimum-Cost Ontology]]
