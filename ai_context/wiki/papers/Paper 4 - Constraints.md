---
type: paper
domain: apf
layer: 4-constraints
created: 2026-04-14
updated: 2026-04-19
sources:
  - Paper_4_Admissibility_Constraints_Field_Content_v2.0.tex
  - Paper_4_Admissibility_Constraints_Field_Content_Supplement_v1.0.tex
---

# Paper 4 - Constraints
## "Admissibility Constraints and the Standard Model Field Content"

**Status:** v2.0 complete (April 19, 2026) — full-scope .tex rebuild with Technical Supplement v1.0. Transitioned from PDF-only (Feb 2026) to 39pp main + 19pp supplement. No longer source-pending. All 115+ theorems in `gauge.py` and `generations.py` now have a prose home.

**Filenames:**
- `Paper_4_Admissibility_Constraints_Field_Content_v2.0.{tex,pdf}` (39pp main paper)
- `Paper_4_Admissibility_Constraints_Field_Content_Supplement_v1.0.{tex,pdf}` (19pp Technical Supplement)

**Core claim:** The Standard Model's gauge group, 45-fermion spectrum, 3-generation structure, weak mixing angle, masses, and mixings follow from finite admissibility alone. Zero free parameters. Combined with Parts I structural constraints (monogamy, area-law, saturation, horizons), the whole Standard Model field content is derived, not chosen.

**Structure (five Parts):**

*Part I — Structural Constraints (§1–§7).* Preserves the Feb 2026 v1.0 content under PLEC framing:
- Non-closure under composition
- Monogamy $T_M$ (biconditional: monogamy $\iff$ finite capacity) — now [P]
- Area-law scaling
- Saturation as feasibility boundary
- Capacity Saturation $T_{4F}$: $N_{\mathrm{gen}} = 3$ from $E(3) = 7 \leq 8 = C_{\mathrm{int}}$ vs. $E(4) = 9 > 8$
- Horizons as persistent saturation
- Structural exclusions

*Part II — Gauge Structure + Field Content (§8–§11).* NEW in v2.0:
- Theorem R (R1/R2/R3) + $L_{\mathrm{gauge\_template\_uniqueness}}$ + $T_{\mathrm{gauge}}$ ($N_c = 3$ by cost minimization)
- Anomaly-freedom + Witten parity + the 4,680-template scan
- $T_{\mathrm{field}}$: 1-of-4,680 unique survivor (Standard Model)
- 45-fermion decomposition ($Q, L, u^c, d^c, e^c$ × 3 generations)
- Capacity partition $C_{\mathrm{total}} = 61 = 45 + 4 + 12$ via $L_{\mathrm{count}}$
- $N_{\mathrm{gen}} = 3$ alternative derivation via capacity-ladder ($T_7$)

*Part III — Mass, Mixing, Weak Mixing Angle (§12–§15).* NEW in v2.0:
- Full $\sin^2\theta_W = 3/13$ derivation ($T_{25a/b}$, $T_{27d}$, $T_{23}$, gate $S_0$)
- Mass structure: capacity-per-dimension + Yukawa bilinear + Schur chain + $m_c$ at 2.6% limit
- Cabibbo $\sin\theta_C = e^{-3/2}$
- CKM matrix structure + CP violation + GCH bridge
- PMNS matrix + neutrino ordering + $m_{\beta\beta}$ prediction + DUNE response

*Part IV — Dark Sector + Predictions + Open Problems (§16–§19).* Expanded from v1.0:
- DM = $C_{\mathrm{ext}}$ geometric identification (not a particle)
- $\Lambda$ = global capacity residual
- Relation to standard quantum constraints (no-cloning, monogamy, Holevo, area laws as pre-quantum)
- 5 explicit falsifiers with observational targets

*Part V — Closure (§20).* Forward refs to Papers 5–7; honest status of the Standard Model derivation claim.

**Key theorems:**
- $T_M$ (Interface monogamy, biconditional) — [P]
- $T_{4F}$ (3 generations) — [P_structural]
- Theorem_R (R1/R2/R3 checks) — [P]
- $L_{\mathrm{gauge\_template\_uniqueness}}$ — [P]
- $T_{\mathrm{gauge}}$ (SU(3)×SU(2)×U(1), $N_c = 3$) — [P]
- $T_{\mathrm{field}}$ (45 fermions, 1-of-4680) — [P]
- $L_{\mathrm{anomaly\_free}}$ — [P]
- $L_{\mathrm{count}}$ ($C_{\mathrm{total}} = 61$) — [P]
- $T_{25a}$, $T_{25b}$ — [P] / [P_structural]
- $T_{\mathrm{sin2theta}}$ ($\sin^2\theta_W = 3/13$) — [P_structural] | gate $S_0$
- $T_{\mathrm{mass\_ratios}}$ — [P_structural]
- $T_{\mathrm{CKM}}$ — [P_structural]
- $T_{\mathrm{PMNS}}$ — [P_structural]

**Paper-specific structural hypotheses (H1–H4), documented in Supplement §1:**
- H1: Generation-channel bridge (GCH) — used for CKM, PMNS, Cabibbo, Schur chain
- H2: Sharing-admissibility (S) — inherited from Paper 2
- H3: Generation cost formula $E(n) = 2n + (n \bmod 2)$ at saturation (alternative $T_7$ ladder form $E(N) = N\varepsilon + N(N-1)\eta/2$ agrees at $N_{\mathrm{gen}} = 3$)
- H4: Gate $S_0$ (interface-schema invariance) — load-bearing for the $3/13$ synthesis

**Implementation:** `gauge.py` (29 checks) + `generations.py` (86 checks). All theorems anchored via `\coderef{check_X}{module.py}`.

**Technical Supplement (v1.0) content:**
- §1 Assumption inventory (PLEC, Regime R, H1–H4)
- §2 $T_M$ biconditional proof (finite $\iff$ monogamy)
- §3 $T_{4F}$ full derivation including $C_{\mathrm{int}} = 8$ justification and CP-coherence correction at $N_{\mathrm{gen}} \geq 4$ (explains why $E(4) = 9$ not $8$)
- §4 1-of-4680 Phase 1 (7 filters, 4680 candidates) + Phase 2 (5 closed-form exclusions)
- §5 $\sin^2\theta_W = 3/13$ gate $S_0$ synthesis (with $4/21$ vs $3/13$ resolution)
- §6 FN ladder + Yukawa bilinear explicit form
- §7 PMNS structural derivation
- §8 Six countermodels (one per PLEC component + one per H-hypothesis)
- §9 Red-team analysis of H1–H4 attack surfaces
- §10 Theorem index + ASCII dependency diagram
- §11 Changelog from Feb 2026 v1.0

**Falsifiers:**
- 4th-generation fermion discovered → invalidates $T_{4F}$, H3
- Direct-detection dark-matter signal → invalidates $C_{\mathrm{ext}}$ geometric identification
- $\sin^2\theta_W$ drift outside $[0.2300, 0.2316]$ at improved precision → invalidates $T_{25b}$, $T_{27d}$, or gate $S_0$
- Non-flat universe ($\Omega_m + \Omega_\Lambda \neq 1$) → invalidates capacity-partition closure
- Anomaly-free fermion template outside the 4,680 scanned → invalidates 1-of-4680 uniqueness

**Open problems:**
- Absolute mass scales $m_t, m_b$ ($\sigma$ derivation required; Paper 7 develops the internalization framework)
- $m_c$ at 2.6% — framework's structural limit at current derivation depth; irreducible error pending deeper Schur-chain corrections
- Neutrino mass hierarchy — DUNE/JUNO will resolve
- DESI-era dynamical dark energy — watched; DR2/DR3 will clarify against the framework's structural $w = -1$ commitment

**Status notes:**
- v1.0 (Feb 2026): 9-page PDF-only. Focus: monogamy, saturation, horizons, sin²θ_W, dark sector.
- v2.0 (April 19, 2026): 39pp main + 19pp supplement. Full-scope rebuild absorbing v1.0 content + gauge derivation + fermion spectrum + mass/mixing previously carried only in codebase.
- Feb 2026 PDF archived as `Old/Paper_4_...vPDFonly_pre-v2.0-rebuild.pdf`.

**Changelog:**
- April 19, 2026 — v1.0 → v2.0. Full-scope rebuild. Five new Parts. Technical Supplement v1.0 created. All 115+ codebase theorems in `gauge.py` + `generations.py` now have a prose home. Released to GitHub (commit `49839a9`).

## See also
- [[Gram Matrix]]
- [[Predictions Catalog]]
- [[Paper 2 - Structure]]
- [[Paper 7 - Action]] (for Type-I seesaw and $\sigma$ derivation)
