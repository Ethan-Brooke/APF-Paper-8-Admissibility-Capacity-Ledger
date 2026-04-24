# Do Not Claim — Paper 8

Anti-hallucination guard. Each row below is a specific negative assertion with a citation to why it is a trap. If you are about to write something that matches a row here, stop.

These rows are harvested from (a) external reviewer responses to Paper 8 v2.0 → v2.5, (b) internal audit findings (the 2026-04-23 smuggling audit), (c) categorical mismatches seen in AI-generated summaries of APF papers. This is not defensive — it is effective steering.

## Core overclaim traps

1. **Do not claim Paper 8 derives $SU(3) \times SU(2) \times U(1)$.** Paper 8 imports that result from Paper 2 (Theorem R + $L_{\rm gauge\_template\_uniqueness}$). *Source: review feedback noted that reviewers repeatedly read the §1.2 irrep table as a derivation.*

2. **Do not claim $I_1$ is a new derivation of the Bekenstein–Hawking entropy.** The cell-count definition $K_{\rm horizon} = A / (4\ell_P^2)$ is a reparametrisation of the standard result; Theorem H.5 Proposition 4.3 establishes equivalence, not novelty. *Source: v2.0 reviewer explicitly asked if this was a new BH calculation; v2.1 Sprint 1 rewrote to clarify it was not.*

3. **Do not claim the Bayes factor in §I'.4 is a full Planck likelihood analysis.** It is a restricted sanity check. It does not marginalise over nuisance, does not use the actual Planck likelihood, does not run MCMC/nested sampling, does not combine with BAO/SN, and does not release chains. Prior-sensitivity table (flat / 0.1-Gaussian / Planck-matched → 0.7 / 1.1 / 2.3) is order-unity and explicitly disclosed. *Source: v2.4 Sprint D reframed §I'.4 as "restricted sanity check, not a full cosmological likelihood analysis".*

4. **Do not claim the repo proves the full 9-paper APF framework.** Paper 8 proves the ACC-ledger consistency architecture and Theorem 1.1. Every other spine result (gauge derivation, field content, T12 partition, $d_{\rm eff} = 102$, Born rule, etc.) is imported from Papers 1–7. *Source: AI summaries sometimes conflate the repo's verification output with a proof of the entire framework.*

5. **Do not claim $d_{\rm eff} = 102$ is derived from scratch in Paper 8.** It is imported from Paper 6's $T_{\rm horizon\_reciprocity}$ + $L_{\rm self\_exclusion}$. Paper 8 uses it. *Source: smuggling-audit finding — the T11 route "27 + 3 + 12" is independent, but it lives in `apf/gravity.py`, not Paper 8.*

6. **Do not claim the $H_0$ tension is resolved here.** Paper 8's $H_0 = 70.03$ km/s/Mpc sits at the midpoint of the Planck–SH0ES tension. The paper states the 7.09$\sigma$ tension as a **falsifier** in Paper 6 §11.4, not as a resolution. *Source: H0DN 2026 framing carried forward from Paper 6.*

7. **Do not treat `[P_structural]` as `[P]`.** $[P_{\rm structural}]$ is conditional on an upstream derivation chain; $[P]$ requires every step to be proved from A1. They are distinct epistemic statuses. *Source: smuggling audit — some numerical-formula checks are `[P_structural]` or `[P_arith]`, not `[P]`.*

## Specific-number traps

8. **Do not claim $C_{\rm vacuum} = 42$ is derived in Paper 8.** It is derived in `apf/gravity.py` T11 via the 27 + 3 + 12 decomposition, and in `apf/cosmology.py` L_equip via $61 - 3 - 16 = 42$. Paper 8 consumes both. *Source: smuggling audit — a reviewer asking "where does 42 come from in Paper 8?" should be answered "imported from T11 + L_equip".*

9. **Do not cite Theorem 1.1 as a proof of $42/61$ without the T12 partition import.** The Maschke argument forces uniqueness of $V_\Lambda$ **given** the irrep-level specification of $V_{\rm local}$. If $V_{\rm local}$ is treated as set-theoretic rather than irrep-content-fixed, the proof flattens. *Source: smuggling audit — this is the strongest attack surface on the paper.*

10. **Do not present the $3 + 16$ subdivision of $V_{\rm local}$ as unquestioned.** "Baryonic = 3 = number of generations" and "Dark = 5 multiplets × 3 gens + 1 Higgs = 16 (including Higgs)" are identifications imported from L_equip. These are contestable — the dark sector inclusion of Higgs is not self-evident. *Source: smuggling audit.*

11. **Do not claim the $42/102^{62}$ coefficient is numerically unique.** An exhaustive scan in `apf/lambda_absolute.py` finds 20 APF-native candidates within 0.01 decades; $19/45$ is numerically closer than $42/102$. The coefficient is privileged **structurally** (two independent bank-forced readings), not numerically. The honest status is `[C]` on the structural-coefficient derivation. *Source: `lambda_absolute.py` coefficient-degeneracy audit.*

12. **Do not read the numerical match $10^{-122.91} \approx 10^{-122.90}$ as evidence that outranks the structural derivation.** The match is a consequence of the structural core, not evidence **for** it. A reviewer who accepts the numerical match but rejects the structural argument has not accepted Paper 8. *Source: Phase 14e.5 retraction — a numerological 12/13 fit was registered and retracted within a single session because it improved the match but not the structural coherence.*

## KMS / Type III$_1$ traps

13. **Do not claim the Type III$_1$ factor classification is novel to APF.** It is a direct consequence of Araki 1968 + Kastler–Haag + Connes 1973 once the algebra, dynamics, and state are specified. What is APF-specific is (i) the UHF algebra indexed by ledger slots, (ii) the vacuum-residual projector structure forced by Theorem 1.1, (iii) identification of partition function and expectation formulas with ledger projections $\pi_A, \pi_Q, \pi_T$. *Source: v2.4 Sprint C added a "Standard imported machinery vs APF-specific content" block at the top of §H.5.*

14. **Do not misstate the assumptions of Theorem 7.14'.** A1 (non-degenerate local Hamiltonian with ≥ 2 distinct eigenvalues), A2 (translation action), A3 (limit state existence). A single-eigenvalue local Hamiltonian would fail Step 4 of the proof and give Type II$_1$ or similar. *Source: v2.4 Sprint B added explicit Assumptions block and "What could go wrong" remark.*

## Minimal working example traps

15. **Do not treat the $K = 3$, $d_{\rm eff} = 4$ toy as a physical system.** It is a pedagogical instance demonstrating the ledger machinery's well-definedness at small scale. Any attempt to "derive physics from the toy" is a category error. *Source: Supplement §R.9 makes this explicit.*

16. **Do not claim the Cobaya YAML stub is runnable as committed.** It is a 20-line reference showing what a community-standard pipeline would look like, not an executable artifact. Running it requires likelihood + MCMC infrastructure Paper 8 does not ship. *Source: v2.4 Sprint D.*

## Scope-extension traps

17. **Do not extend scope from "consistency of existing PLEC-chain derivations" to "microscopic derivation".** §N.6 explicitly owns the "consistency reparametrisation" characterisation. APF is a meta-theoretical consistency framework complementary to microscopic theories (LQG, string, SUSY-localisation), not a replacement. *Source: v2.2 Sprint 10.*

18. **Do not claim the repo contains the full APF codebase.** The bundled `apf/` subset is a task-scoped extract for Paper 8 verification. Full codebase v6.9 contains 420 bank-registered theorems across 34 modules; this repo contains the subset needed to verify Paper 8's claims. *Source: `START_HERE.md` §0 already makes this explicit.*

19. **Do not treat `§M.7` modification-catalogue rejections as proofs.** They are observational exclusions (what Planck 2018 + FIRAS + SH0ES rule out), not structural impossibilities. *Source: v2.2 Sprint 1.*

20. **Do not claim `K_horizon` is observer-invariant.** Proposition 16.6 establishes $K_{\rm horizon}$ is observer-covariant with transformation law $K_{\rm horizon}(o') = A(o') / (4\ell_P^2)$; only $K_{\rm SM}$ and $d_{\rm eff}^{\rm SM}$ are observer-invariant (Propositions 16.5 and 16.7). *Source: v2.2 Sprint 7.*

---

*20 rows. Review this file before writing any summary, email, or reviewer response that mentions Paper 8. If you find a predictable AI overclaim not covered here, add a row.*
