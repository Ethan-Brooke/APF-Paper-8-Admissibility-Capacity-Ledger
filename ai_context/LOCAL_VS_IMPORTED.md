# Local vs Imported — Paper 8

This file partitions every structurally meaningful result in Paper 8 into one of five categories. Its purpose: prevent the most common AI failure mode, which is presenting imported APF-spine results as if Paper 8 derives them from scratch.

If you are a cold AI onboarding this repo, read this file before attempting to summarise the paper. If you find yourself about to write "Paper 8 derives the SM gauge group" or "Paper 8 proves $d_{\rm eff} = 102$", stop and check the table below.

## Locally proved in Paper 8

- Formal ACC record $(K, \mathrm{ACC})$ as a two-scalar object with explicit algebraic structure.
- Six regime projections $\pi_T, \pi_G, \pi_Q, \pi_F, \pi_C, \pi_A$ and their well-definedness on the ACC record.
- Three canonical interface factories `acc_SM`, `acc_horizon`, `acc_quantum`.
- The four-identity architecture $I_1, I_2, I_3, I_4$ as a stateable, falsifiable collection of consistency theorems.
- **Theorem 1.1 (gauge-cosmological bridge) in Supplement §1.4:** existence and uniqueness of a 42-dim $G_{\rm SM}$-invariant $V_\Lambda \subset V_{61}$, via Maschke complete reducibility.
- Integer-level closure $3 + 16 + 42 = 61$ as a corollary of Theorem 1.1 (demoted to §D.2 in v2.5).
- **Theorem 5.6 (SM ACC record uniqueness)** under four stated upstream identifications (U1)–(U4).
- Theorem 7.14' (Type III$_1$ factor at infinite volume) as a composition of the listed standard imports; APF-specific content is only the choice of algebra and the vacuum-residual projector structure.
- Proposition 3.4 (two-tier consistency forcing) reducing five integer commitments to two at the SM interface.
- The minimal working example at $K=3$, $d_{\rm eff}=4$ (Supplement §R): all six projections evaluated in closed form, one bridge identity verified by hand, one operator expectation verified, 10-line numpy reproduction.
- Prior-sensitivity table for the restricted Bayes-factor sanity check in §I'.4.
- Meta-theoretical reframing in §N.6 (APF as consistency-reparametrisation framework).

## Imported from earlier APF papers

- Gauge group $SU(3) \times SU(2) \times U(1)$ derivation — **Paper 2**, via Theorem R + $L_{\rm gauge\_template\_uniqueness}$.
- SM field-content count $12 + 4 + 45 = 61$ — **Paper 2 / Paper 4**, via $L_{\rm count}$ + the 1-of-4680 fermion-content filter.
- T12 partition constraint $V_{61} = V_{\rm local} \oplus V_\Lambda$ with $\dim V_{\rm local} = 19$ — **Paper 6 Supplement §T12**.
- $d_{\rm eff} = 102$ via $T_{\rm horizon\_reciprocity}$ + $L_{\rm self\_exclusion}$ — **Paper 6 Supplement**.
- $C_{\rm vacuum} = 42$ with structural decomposition $27 + 3 + 12$ — **Paper 6 Supplement §T11**, independent of the matter-sector subtraction route.
- Cosmological-partition motivation for the $\Omega_\Lambda = 42/61$ reading — **Paper 6 Supplement §T11** and **Paper 3 Supplement §7** (T_ACC formalization).
- PLEC four-component framework (A1 + MD + A2 + BW) — **Paper 1 v4.0-PLEC**.
- Matter-sector equipartition $n_{\rm baryon} = 3$, $n_{\rm dark} = 16$, $n_{\rm vacuum} = 42$ — **Paper 6 / apf/cosmology.py L_equip**. Paper 8 consumes these; it does not derive the identification.

## Standard mathematical imports

- **Maschke complete reducibility** for compact Lie groups. Cited: Hall 2015. Used in Theorem 1.1 proof.
- Finite-dimensional von Neumann entropy closure: $S_{\rm vN}(\rho_{\max\text{-mixed}}) = \ln \dim \mathcal{H}$. Textbook.
- High-temperature partition-function limit $\ln Z(\beta) \sim K \ln d_{\rm eff}$ as $\beta \to 0$. Textbook.
- $C^*$-inductive-limit machinery (Bratteli–Robinson) where used in §H.5.
- Araki 1968 extremal-KMS factor property; Kastler–Haag finite-volume factor property; Connes 1973 Type III$_1$ classification. Composed in Theorem 7.14' proof.
- Gleason's theorem (cited in main paper §5 via Paper 5 derivation of the Born rule).
- Lieb–Ruskai strong subadditivity (via Paper 3 entropy backbone).
- Tikz-cd for the bridge-theorem commutative diagram (Supp §E.3) — visual aid, not a proof step.

## Empirical comparisons only (not derivations)

- **Planck 2018** $\Omega_\Lambda$ posterior $0.6847 \pm 0.0073$ — used in the §I'.4 restricted sanity check.
- **Planck 2018** $3 \times 3$ covariance matrix over $(\Omega_b, \Omega_c, \Omega_\Lambda)$ with $\rho_{c\Lambda} \approx -0.95$ — used for the $\chi^2 = 1.18$ multivariate sanity check.
- **SH0ES 2022** $H_0 = 73.04 \pm 1.04$ km/s/Mpc — used only in the H_0-tension framing. The 7.09$\sigma$ tension is stated as a **falsifier**, not resolved by the paper.
- **FIRAS** $T_{\rm CMB}$ measurement — used in the temperature-formula agreement statement.
- **Riess 2022** $H_0$, used in the brittle-consequences table in §I'.
- Ong et al. 2511.10631, 2511.04661 — recent Bayesian re-analyses cited as context for §I'.4's prior-sensitivity discussion.

## Not claimed by this paper (common overclaim risks)

- **A full cosmological likelihood analysis.** §I'.4 is explicitly a restricted sanity check; Cobaya stub is 20 lines and not runnable as committed.
- **Microscopic horizon counting.** Paper 8 does not compete with LQG puncture counting, string D-brane state counting, or SUSY-localisation horizon entropy calculations. §N.1.1 positions APF as complementary, not competing.
- **A new derivation of Bekenstein–Hawking entropy.** The cell-count definition $K_{\rm horizon} = A / (4\ell_P^2)$ is a reparametrisation of the standard result, not a new calculation.
- **Independent derivation of all APF spine results.** Paper 8 imports the gauge group, field content, T12 partition, and $d_{\rm eff} = 102$ from Papers 1/2/4/6. It formalises their joint consistency at the ledger level; it does not re-derive any of them.
- **Resolution of the $H_0$ tension.** Paper 6 states the tension as a falsifier; Paper 8's $H_0 = 70.03$ km/s/Mpc prediction sits at the tension midpoint but does not explain either endpoint.
- **An operational procedure for measuring $\varepsilon^*$ in a lab.** §A.4.1 gives calibration protocols at three interfaces; an absolute $\varepsilon^*$ measurement is explicitly flagged as open in §A.4.2.
- **Structural uniqueness of the $\rho_\Lambda / M_{\rm Pl}^4 = 42/102^{62}$ formula against all numerical competitors.** `apf/lambda_absolute.py` exhaustively scans and finds 20 APF-native candidates within 0.01 decades; $19/45$ is numerically closer. The chosen coefficient is privileged **structurally** (two independent bank-forced readings via $L_{\rm self\_exclusion}$ and $T_{12}$ + $T_{\rm interface\_sector\_bridge}$), not numerically. The paper is honest about this: `[C]` status on the structural-coefficient derivation.
- **That $[P_{\rm structural}]$ has the same epistemic status as $[P]$.** It does not. Structural derivations are conditional on the upstream derivation chain; $[P]$ requires every step to be proved from A1.

## How to use this file during onboarding

If you are preparing to answer any of the following, consult this file first:

- "What does Paper 8 derive?" → **Locally proved** section.
- "What does Paper 8 assume?" → **Imported** and **Standard mathematical imports** sections.
- "Is this a full likelihood analysis?" → See **Not claimed**.
- "Does the paper resolve the Hubble tension?" → See **Not claimed**.
- "Where does $d_{\rm eff} = 102$ come from in this paper?" → **Imported**, not locally proved.

If a claim appears in the supplement or main paper and is *not* traceable to one of these five categories, flag it as a categorisation gap.
