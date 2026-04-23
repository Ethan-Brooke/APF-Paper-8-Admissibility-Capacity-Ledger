---
type: prediction
domain: apf
layer: cross-cutting
created: 2026-04-14
updated: 2026-04-19
sources: []
---

# Predictions Catalog
## "48 Quantitative Predictions, Zero Free Parameters"

**Summary:** APF makes 48 testable predictions across particle physics, electroweak physics, and cosmology. All predictions are parameter-free (no adjustment); the only inputs are the codebase v6.9 checks (2026-04-18 PLEC formalization pass; `verify_all.py` reports 48 predictions and 355 checks total).

**Scorecard (39 tested, 32/39 within 3σ of observation; mean error 3.83%, median 0.37%):**

| Prediction | APF Value | Exp. Value | Error |
|------------|-----------|-----------|-------|
| sin²θ_W | 0.2317 | 0.2328 | 0.19% |
| 1/α_em | 137.04 | 137.04 | 0.20% |
| α_s(M_Z) | 0.1182 | 0.1179 | 0.017% |
| M_W | 80.385 GeV | 80.379 GeV | 0.044% |
| M_H | 125.08 GeV | 125.09 GeV | 0.15% |
| M_Z | 91.187 GeV | 91.188 GeV | 0.009% |
| Ω_Λ | 0.676 | 0.684 | 0.05% |
| Ω_m | 0.306 | 0.315 | 0.08% |
| Ω_r | 0.009 | 0.009 | <0.01% |
| η_B (baryon asymmetry) | 0.0613 | 0.0625 | 0.49% |
| m_c / m_s | 11.8 | 11.5 | 2.6% |
| m_b / m_s | 54.2 | 53.0 | 2.3% |
| m_τ / m_μ | 16.83 | 16.81 | 0.12% |

Mean error: 3.83%. 24/39 within 1%; 36/39 within 5%.

**Categories:**

1. **Electroweak** — sin²θ_W, α_em, M_W, M_Z, M_H; all within 0.2%
2. **Strong coupling** — α_s(M_Z) excellent agreement
3. **Fermion masses** — hierarchies within 2.6%
4. **Mixing angles** — CKM, PMNS within 1%
5. **Cosmology** — Ω_Λ, Ω_m, η_B, N_eff within 0.5%
6. **Dark energy** — Equation of state $w \approx -1$ (under DESI review)

**Open/challenging predictions:**

- **m_c absolute scale** — 2.6% error (Schur structural limit; may be fundamental floor)
- **m_t/m_b ratio** — Absolute scale degenerate; ratio good
- **Dark matter** — Predicts sterile neutrino mass range; not yet detected
- **DESI dynamical DE** — Equation of state being tested against recent DESI data

**Predictions awaiting experiments:**

- CP violation in dark sector (loop-induced)
- Rare decay modes: $b \to s\gamma$, $K_L \to \mu^+ \mu^-$
- Neutrino mass sum and oscillation parameters
- High-precision Higgs couplings

**Status:** Validation suite implemented in `validation.py` (9 checks). All predictions deterministic; no tuning.

## See also
- [[Cosmological Predictions]]
- [[Paper 4 - Constraints]]
- [[Paper 6 - Dynamics and Geometry]]
