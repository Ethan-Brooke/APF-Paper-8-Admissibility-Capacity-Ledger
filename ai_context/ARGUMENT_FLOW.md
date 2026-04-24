# Argument Flow — Paper 8 in One Page

A distilled walk-through of Paper 8's argument end to end. Not a DAG, not a README, not the supplement's reader map. This is the structural spine: read it in order, and you have the paper.

**If you are reading only one file in `ai_context/`, read this one.**

## A. Ledger architecture

A1 (finite capacity, imported from Paper 1) is the single axiom. PLEC (Principle of Least Enforcement Cost) gives four structurally necessary components: A1 + MD + A2 + BW.

The **ACC record** is a two-scalar object $(K, \mathrm{ACC})$ where:

- $K$ is the structural capacity — integer count of admissibility channels at an interface.
- $\mathrm{ACC} = K \ln d_{\mathrm{eff}}$ is the log-count scalar; $d_{\mathrm{eff}}$ is the per-slot effective dimension.

Six **regime projections** read the ACC record in six physical modes: $\pi_T$ (thermo), $\pi_G$ (gravitational), $\pi_Q$ (quantum), $\pi_F$ (gauge/field), $\pi_C$ (cosmological), $\pi_A$ (action).

Three **canonical interface factories**: `acc_SM(K=61, d_eff=102)` for the Standard Model interface, `acc_horizon(A)` for a horizon of area $A$, `acc_quantum(d)` for a quantum interface of dimension $d$.

## B. Four identities

The ACC machinery has to be consistent across the six regimes. Four consistency identities enforce this:

- $I_1$ — **horizon convention identity.** $K_{\rm horizon}$ read in nats (Convention A) and in bits (Convention B) agree at the Bekenstein–Hawking-matching interface. *Status: definitional.*
- $I_2$ — **gauge-cosmological bridge.** The integer $K$ that $\pi_F$ reads at the SM interface equals the denominator that $\pi_C$ uses to form $\{\Omega_b, \Omega_c, \Omega_\Lambda\}$. Both are 61. *Status: nontrivial (this is the paper's structural core).*
- $I_3$ — **thermo-quantum closure.** On a max-mixed state, $S_{\rm vN} = \ln \dim \mathcal{H} = \mathrm{ACC}$. *Status: standard finite-dim result, locally proved.*
- $I_4$ — **action-thermo closure.** $\lim_{\beta \to 0} \ln Z(\beta) = K \ln d_{\rm eff}$. *Status: standard partition-function limit, locally proved.*

**Only $I_2$ is load-bearing.** $I_1, I_3, I_4$ are modest consistency identities. This hierarchy is explicit in Supp §1.6 and in the paper's abstract.

## C. Nontrivial core (Supplement §1 — "Formal Kernel of the Paper")

Front-loaded immediately after the reader map. Six subsections, one theorem, one proof.

1. **Define** $V_{61}$ — a 61-dimensional complex vector space with explicit ordered basis $\{e_i^{(\rm gauge)}\}_{i=1}^{12} \cup \{e_j^{(\rm Higgs)}\}_{j=1}^{4} \cup \{e_k^{(\rm fermion)}\}_{k=1}^{45}$. Dimension $12 + 4 + 45 = 61$.

2. **Define** $G_{\rm SM} = SU(3) \times SU(2) \times U(1)$ acting on $V_{61}$ via representation $\rho_{\rm SM}$. Full irrep decomposition given: adjoint $8 + 3 + 1 = 12$ for gauge; complex doublet $Y = 1/2$ giving 4 for Higgs; $Q_L + u_R + d_R + L_L + e_R$ per generation = $15 \times 3 = 45$ for fermion.

3. **Define** the T12 partition constraint $V_{61} = V_{\rm local} \oplus V_\Lambda$ with $\dim V_{\rm local} = 19$, $V_{\rm local} = V_b \oplus V_c$ with $\dim V_b = 3$, $\dim V_c = 16$.

4. **Prove** (Theorem 1.1): there exists a unique $G_{\rm SM}$-invariant 42-dim subspace $V_\Lambda$ with $V_{61} = V_{\rm local} \oplus V_\Lambda$. Via **Maschke complete reducibility** for compact Lie groups (Hall 2015 cited as the only external dependency). Four proof steps, no APF-specific machinery as proof step.

5. **Consequence:** $\Omega_\Lambda = \dim V_\Lambda / \dim V_{61} = 42/61$ under the Fractional Reading Equivalence.

**Honest uncertainty at this step.** The Maschke argument forces uniqueness of $V_\Lambda$ **given the irrep specification of $V_{\rm local}$**. If $V_{\rm local}$ is specified only by "everything that is not $V_\Lambda$", the proof is circular. The paper's resolution: Theorem 1.1 imports $V_{\rm local}$'s irrep content from Paper 4 / Paper 6, not from Paper 8. See `LOCAL_VS_IMPORTED.md` for the boundary.

**Executable status (Phase 22.2, in progress).** Paper proof is complete. Standalone code witness `check_T_FormalKernel_VLambda_uniqueness` is landing Phase 22.2.a (next session) — constructs $V_{61}$ with explicit $G_{\rm SM}$ matrix reps, enumerates $V_{\rm local}$ by irrep content, computes the invariant complement, and verifies dim = 42 plus uniqueness via an adversarial dimension-only counterexample. Currently the downstream witness `check_T_interface_sector_bridge` certifies the $V_{\rm global} \cong V_\Lambda$ identification at the DAG level.

## D. Downstream predictions

Using $K = 61$, $d_{\rm eff} = 102$, and the $3 + 16 + 42$ partition:

- $\Omega_b = 3/61$, $\Omega_c = 16/61$, $\Omega_m = 19/61$, $\Omega_\Lambda = 42/61$. Planck 2018 agrees to 0.5$\sigma$ (on $\Omega_\Lambda$) and better on matter-sector ratios.
- $\rho_\Lambda / M_{\rm Pl}^4 = 42 / 102^{62}$. Numerical value $10^{-122.91}$; observed $10^{-122.90}$; agreement 3% (= 0.012 decades). **Structural coefficient is `[C]`** — an exhaustive scan finds near-neighbour competitors (notably $19/45$), so the coefficient is privileged structurally (two independent bank-forced readings), not numerically.
- $H_0 = 70.03$ km/s/Mpc via GR $\rho_{\rm crit}$ inversion. Sits at the midpoint of the Planck–SH0ES tension; the paper states the 7.09$\sigma$ tension as a **falsifier**, not a resolution.
- $(T_{\rm CMB} / M_{\rm Pl})^4 = 48 / 102^{64}$. Matches FIRAS at 0.33%.
- Two-tier-structure brittleness: §I' sensitivity theorem bounds how much any of these predictions can shift if upstream integer counts change by 1.

All empirical comparisons are explicitly framed as consequences of the structural core, not evidence **for** it. The structural core stands or falls on Theorem 1.1, not on $H_0$.

## E. Minimal working example (Supplement §R)

A pedagogical toy interface at $K = 3$, $d_{\rm eff} = 4$, $C_{\rm vacuum} = 1$. Fully self-contained:

- $V_3 = \mathbb{C}^3$, trivial gauge group, residual partition $(1, 1, 1)$.
- $H_{\rm micro} = (\mathbb{C}^4)^{\otimes 3} = \mathbb{C}^{64}$, algebra $M_{64}(\mathbb{C})$.
- All six projections in closed form: $\pi_F = 3$, $\pi_Q = 64$, $\pi_T = 3 \ln 4$, $\pi_C = (1/3, 1/3, 1/3)$, etc.
- One bridge identity verified by hand: $\pi_F = 3 =$ denominator of $\pi_C$.
- One operator expectation verified: $\langle P_{\rm vac} \rangle = 16/64 = 1/4 = C_{\rm vac}/d_{\rm eff}$.
- 10-line numpy reproduction: `minimal_working_example/toy_interface_numpy.py`.

The MWE teaches the **shape** of the argument without requiring APF-specific vocabulary. Running it and comparing to the committed expected output is a 30-second audit of the ledger machinery.

---

**If you just read this file, you understand the paper.** Next: read `LOCAL_VS_IMPORTED.md` to learn which steps are Paper 8's own proofs and which are imported; read `DO_NOT_CLAIM.md` before writing about the paper; read `CLAIMS_LEDGER.md` for the row-by-row attack surface.

*One page, no exceptions. If this file grows beyond a single print page, it has lost its purpose.*
