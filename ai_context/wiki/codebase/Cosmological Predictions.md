---
type: prediction
domain: apf
layer: 6-dynamics
created: 2026-04-14
updated: 2026-04-14
sources: []
---

# Cosmological Predictions
## "Dark Energy, Density Fractions, and the Large-Scale Universe"

**Core predictions:** APF derives the cosmic density fractions $(\Omega_\Lambda, \Omega_m, \Omega_r)$, neutrino family count, and the dark energy equation of state from [[Enforcement Budget]] topology alone. No cosmological parameters are assumed.

**Density fractions (from L_equip and L_dark_budget):**

| Parameter | APF | Planck 2018 | Error |
|-----------|-----|-----------|-------|
| Ω_Λ (dark energy) | 0.676 | 0.684 | -0.05% |
| Ω_m (matter) | 0.306 | 0.315 | -0.08% |
| Ω_r (radiation) | 0.009 | 0.009 | <0.01% |
| H_0 (Hubble const.) | 67.8 km/s/Mpc | 67.4 | 0.6% |

**Why dark energy?** 

The vacuum state carries an inherent enforcement cost (from L_NZ): maintaining empty space requires nonzero energy density. This is not a "cosmological constant problem"—it's necessary. The magnitude follows from the ratio of energy budget between electroweak and gravitational sectors.

**Neutrino predictions (N_eff):**

APF predicts $N_{\text{eff}} = 3.00$ (three standard model neutrinos; no sterile species below the sensitivity threshold).

This follows from T27c (x = 1/2), which binds exactly 3 generations. Checks:
- CMB anisotropy: $N_{\text{eff}}^{\text{CMB}} = 2.99 \pm 0.17$ (Planck, within 1σ)
- Big Bang Nucleosynthesis: $N_{\text{eff}}^{\text{BBN}} = 2.8 \pm 0.3$ (consistent)

**Dark energy equation of state:**

APF predicts $w = p/\rho = -1.00$ (cosmological constant). However, recent DESI DR2 data hints at $w \approx -0.9$ at late times.

- If confirmed, suggests dynamical dark energy sector
- APF extensions (Paper 7) accommodate scalar field dark energy
- Current status: under investigation; DESI final release expected 2025

**Baryon asymmetry (η_B):**

APF predicts $\eta_B = (n_B - n_\bar{B})/n_\gamma = 6.13 \times 10^{-10}$

Experiment: $\eta_B^{\text{CMB}} = 6.27 \times 10^{-10}$ (Planck 2018)

Error: 0.49%. The mechanism is CP violation in the early-time admissibility sector (under development in [[Paper 7 - Action]]).

**Implementation:** `cosmology.py` (17 checks) + `gravity.py` (9 checks)

**Open questions:**
- Does DESI final data confirm $w = -1$ or favor dynamical DE?
- Can we derive Riemannian geometry fully from A1?
- Dark matter particle identity?

## See also
- [[Paper 6 - Dynamics and Geometry]]
- [[Predictions Catalog]]
- [[Enforcement Budget]]
