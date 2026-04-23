# FRAMEWORK_OVERVIEW.md — APF in 5 minutes

This file is the framework-wide context that one paper's companion repo cannot supply on its own. Read after `AGENTS.md`.

---

## The axiom

**A1 (Finite Information Capacity):** Any interface of physical distinction has finite information capacity.

Formally: for any interface $\Gamma$ with resource structure $(\rho, \Gamma)$, there exists a finite $C(\rho, \Gamma)$ such that the sum of enforcement costs of all distinctions maintained at $\Gamma$ is bounded by $C$:

$$
\sum_{d \in G} \varepsilon(d) \leq C(\rho, \Gamma).
$$

Everything in APF is derived from A1 plus three structural companions:

## The four-component PLEC (Principle of Least Enforcement Cost)

A1 alone bounds the admissible set from above. Three additional components make the structure non-trivial:

| Component | Role | Content |
|---|---|---|
| **A1** (finite capacity) | upper bound | $\sum \varepsilon \leq C$ |
| **MD** (minimum distinction cost) | lower bound | $\varepsilon(d) \geq \mu^* > 0$ for all $d$ |
| **A2** (argmin selection) | realization | physical structure is least-cost member of admissible class |
| **BW** (non-degeneracy) | spectrum | cost spectrum is non-degenerate |

Together these four components are a single variational principle — PLEC — that selects realized physical structure from the admissible class.

The four components are independent: each has a countermodel where the other three hold but it fails. See `apf/core.py:check_L_structural_independence`.

## The admissibility structure

The object APF actually reasons about:

$$
\mathcal{A}(\rho, \Gamma) = \left\{ G \mid \varepsilon: G \to \mathbb{R}_{>0},\; \inf_{d \in G}\varepsilon(d) \geq \mu^*,\; \sum_{d \in G}\varepsilon(d) \leq C(\rho, \Gamma),\; \varepsilon\text{ non-degenerate}\right\}.
$$

Bounded by MD from below, A1 from above, BW for spectrum discipline. PLEC (via A2) selects the realized element.

## The 61-partition

The single most important quantitative output of the framework is the type-count partition:

$$
\boxed{C_{\mathrm{total}} = 61 = \underbrace{42}_{\text{vacuum}} + \underbrace{19}_{\text{matter}} = 42 + \underbrace{3}_{\text{baryon}} + \underbrace{16}_{\text{cold dark}}}
$$

The numbers are not tuned. They are forced by:

- **Anomaly-free field content** (L_anomaly_free): the Standard Model gauge group $SU(3) \times SU(2) \times U(1)$ with chiral fermions admits exactly 61 type-classifications satisfying the 7 independent anomaly conditions.
- **Gauge-addressability** ($Q_1$, from T3): splits 61 = 42 + 19 (vacuum types that don't route through non-trivial gauge channels vs. matter types that do).
- **Confinement** ($Q_2$, from T_confinement): splits 19 = 3 + 16 (gauge-addressable types carrying $SU(3)$ color vs. colorless).

The partition is **topological** (type-counting), not dynamical. `check_L_saturation_partition` proves it is saturation-independent: as capacity fills or drains, the ratios stay fixed.

## The derivation chain (one-page sketch)

```
A1 (finite capacity)
 ├─→ L_nc       → non-abelian carrier required
 │                → Theorem_R(R1) → SU(N_c), N_c ≥ 3
 │
 ├─→ L_irr      → chiral carrier (locally unrecoverable correlations)
 │                → Theorem_R(R2) → SU(2) [UNIQUE 2D pseudoreal]
 │
 ├─→ L_col      → minimality
 │                → Theorem_R(R3) → U(1), N_c = 3
 │
 │   ├─→ L_gauge_template_uniqueness: the template is FORCED
 │   ├─→ T_gauge: N_c = 3 selected by cost optimization
 │   ├─→ T_field: 45 chiral fermions (1 of 4680 anomaly-free sets survives)
 │   ├─→ T4F: N_gen = 3 exactly (capacity saturation)
 │   └─→ L_count: C_total = 45 + 4 + 12 = 61 [rigid integer]
 │                → 25 Layer II predictions (masses, mixings, cosmology)
 │
 └─→ L_Cauchy_uniqueness: F(d) = d → γ = 17/4 → sin²θ_W = 3/13
```

Horizontal lines above are logical necessity ("forces"); downward branches are dependent derivations.

## The 19 bank modules (v6.9)

| Module | Role |
|---|---|
| `core` | Axiom A1 / PLEC / Hilbert space / Born rule / CPTP. Foundational layer. |
| `gauge` | Gauge structure, non-closure, uniqueness, record-locking (T9 / L3-µ), particle content. |
| `generations` | Fermion spectrum, 3 generations, masses, mixing, texture, seesaw. |
| `spacetime` | Metric structure, $d = 4$, Coleman–Mandula closure, FBC geometry. |
| `gravity` | Gravitational sector, horizon entropy, T7B (metric), T9_grav (Einstein), T10 (κ scaling), A9 closure. |
| `cosmology` | Cosmological constant (T11), dark-energy equation of state, Planck 4-param match, DESI response. |
| `plec` | **v6.9 NEW.** Regime R + five-type regime-exit taxonomy. |
| `validation` | Prediction validation vs. experiments. |
| `supplements` | Gleason internalization, structural lemmas, computational aids. |
| `majorana` | Neutrino physics, Majorana constraints. |
| `internalization` | Spectral action internalization (Connes / NCG). |
| `internalization_geo` | Geometry from internalization (Kolmogorov, Lipschitz, Lovelock, HKM, Malament). |
| `extensions` | BSM physics, dark-matter mechanisms, HKM causal geometry. |
| `red_team` | Adversarial tests for false positives and edge cases. |
| `session_v63c` | v6.3c incremental: neutrino hierarchy, Yukawa spectral. |
| `session_qg` | Quantum gravity / Planck-scale scoping. |
| `session_nnlo` | NNLO corrections. |
| `session_delta_pmns` | PMNS phase derivation progress. |
| `session_cosmo_update` | Cosmological observables refresh. |

**Plus `apf/standalone/`**: four additional `verify_all`-tracked files (L_Cauchy_uniqueness, L_CKM_resolution_limit, phase1_seesaw_closure, phase5_theorem_R_audit) that are not bank-registered but are run by `verify_all.py`.

Totals: **342 bank-registered theorems · 355 `verify_all` checks · 48 quantitative predictions · 0 free parameters.**

## Quantitative predictions — the headline results

### Cosmology (4-parameter match to Planck 2018)

| Quantity | APF | Planck 2018 | Δ |
|---|---|---|---|
| $\Omega_\Lambda$ | $42/61 = 0.6885$ | $0.6848 \pm 0.007$ | $+0.54\%$ |
| $\Omega_m$ | $19/61 = 0.3115$ | $0.3152 \pm 0.007$ | $-1.17\%$ |
| $\Omega_b$ | $3/61 = 0.0492$ | $0.0493 \pm 0.0005$ | $-0.25\%$ |
| $\Omega_{\mathrm{cdm}}$ | $16/61 = 0.2623$ | $0.2645 \pm 0.006$ | $-0.82\%$ |
| $H_0$ | $67.76$ | $67.36 \pm 0.54$ | within $1\sigma$ |

### Standard Model parameters (representative)

- **$\sin^2\theta_W = 3/13$** (from L_Cauchy_uniqueness + $\gamma = 17/4$).
- **Fermion mass ratios** within structural hierarchy. Most within 5% of observation.
- **$N_{\mathrm{gen}} = 3$** exact.
- **CKM and PMNS** with structural constraints matching observed magnitudes; phases in progress.

### Prediction scorecard

39 tested quantitative predictions. 32/39 within 3σ. Mean error 3.83%. Median error 0.37%.

## Open problems — genuinely unsolved

Don't try to solve these casually. Each has its own research program.

- **H0 tension** (7.09σ vs H0DN 2026). APF predicts 67.76 km/s/Mpc; sharp falsifier candidate. See `OPEN_PROBLEMS.md`.
- **m_t, m_b absolute scales.** Structural ratios match but absolute scale requires σ derivation.
- **m_c at 2.6%.** Schur structural limit; irreducible within current framework.
- **Dark matter particle identity.** T12 derives existence + gravitational properties + no gauge charge; no particle species ID. See Paper 6 §13.
- **PLEC existence/uniqueness theorem.** The four-component variational principle is stated but not fully formalized as an argmin existence + uniqueness result. Paper 6 App. C Project C.1.
- **Dynamical DE.** DESI DR2 prefers $w_0 \approx -0.75$; APF requires $w = -1$ exactly. 1.5–2σ tension; falsifier at DR3 ≥5σ.

## Scope — what APF IS and what it ISN'T

**APF IS:** a derivation of the *structure* of the Standard Model and general relativity from A1 + PLEC. It predicts gauge group, fermion content, generation count, capacity partition, Hilbert-space structure, Born rule, Lorentzian metric, Einstein equations, cosmological constant, and four-parameter Planck cosmology — all with zero free parameters.

**APF IS NOT:** a theory of quantum gravity (no explicit unification of GR with gauge sector at Planck scale), a numerical parameter-fitter (absolute mass scales for top/bottom quarks not yet derived), a theory of dark matter identity (gravitational role yes; particle species no), or a replacement for quantum field theory at the operational level (QFT emerges within APF's coherent regime).

---

*This file is bundled with every APF paper release. It is authoritative at release time. For the latest canonical state, see the APF codebase (`__APF Library/Codebase/APF_Codebase_vX.X/`) and the framework's Zenodo deposits.*
