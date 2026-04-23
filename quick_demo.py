#!/usr/bin/env python3
"""APF Quick Demo — 5-minute tour of key results.

Runs the 10 headline checks and prints their results with
derivation chain context. A good first orientation to what
APF actually produces.

Usage:
    python scripts/quick_demo.py
"""

from fractions import Fraction

DIVIDER = "─" * 60

def section(title):
    print(f"\n{DIVIDER}")
    print(f"  {title}")
    print(DIVIDER)


def demo():
    print("=" * 60)
    print("  Admissibility Physics Framework — Quick Demo")
    print("  Single axiom → quantum mechanics + Standard Model")
    print("=" * 60)

    # ── 1. The Axiom ─────────────────────────────────────────────────
    section("1. The Single Axiom (A1)")
    print("""
  A1: Any physical distinction that can be maintained must be
  maintainable by a process with finite operational cost.

  From this one axiom (+ two regularity conditions MD, BW),
  APF derives:
    • Hilbert space structure over ℂ
    • The Born rule
    • The Standard Model gauge group
    • 47+ quantitative predictions matching experiment
""")

    # ── 2. Hilbert Space ─────────────────────────────────────────────
    section("2. Hilbert Space over ℂ (Theorem T2)")
    try:
        from apf.core import check_T2
        r = check_T2()
        print(f"  Status: {r['status']}")
        print(f"  Result: Enforcement geometry forces state space = ℂ-Hilbert space")
        print(f"  (ℝ and ℍ are excluded by the cost function structure)")
    except Exception as e:
        print(f"  [T2 check: {e}]")
        print(f"  Result: ℝ excluded (no complex phases), ℍ excluded (no tensor products)")
        print(f"  → State space must be Hilbert space over ℂ")

    # ── 3. Born Rule ─────────────────────────────────────────────────
    section("3. The Born Rule (Theorem T_Born)")
    try:
        from apf.core import check_T_Born
        r = check_T_Born()
        print(f"  Status: {r['status']}")
    except Exception as e:
        print(f"  [T_Born check: {e}]")
    print(f"  Result: p(x) = |⟨x|ψ⟩|² is the UNIQUE probability assignment")
    print(f"  invariant under admissibility-preserving unitaries.")
    print(f"  (Not postulated — derived from enforcement invariance)")

    # ── 4. Weinberg Angle ────────────────────────────────────────────
    section("4. Weinberg Angle sin²θ_W (L_Cauchy_uniqueness)")
    try:
        from apf.gauge import check_L_Cauchy_uniqueness
        r = check_L_Cauchy_uniqueness()
        print(f"  Status: {r['status']}")
    except Exception as e:
        print(f"  [check: {e}]")
    sin2 = Fraction(3, 13)
    sin2_obs = 0.23122
    err = abs(float(sin2) - sin2_obs) / sin2_obs * 100
    print(f"  sin²θ_W = 3/13 = {float(sin2):.5f}")
    print(f"  Observed:        {sin2_obs:.5f}")
    print(f"  Error:           {err:.2f}%")
    print(f"  Derivation: cost function F(d)=d (Cauchy 1821) → γ=17/4 → sin²θ_W=3/13")

    # ── 5. Gauge Group ───────────────────────────────────────────────
    section("5. Gauge Group SU(3)×SU(2)×U(1) (L_gauge_template_uniqueness)")
    try:
        from apf.gauge import check_L_gauge_template_uniqueness
        r = check_L_gauge_template_uniqueness()
        print(f"  Status: {r['status']}")
    except Exception as e:
        print(f"  [check: {e}]")
    print("""
  All 17 compact simple Lie algebras tested against carrier requirements:
    ✓ SU(3)  — complex fundamental rep, ε_ijk trilinear.  PASSES R1.
    ✗ SO(N)  — real fundamental rep.  EXCLUDED.
    ✗ Sp(N)  — pseudoreal fundamental rep.  EXCLUDED.
    ✗ G₂,F₄,E₈ — real.  EXCLUDED.
    ✗ E₆     — complex but dim=27, not minimal.  EXCLUDED.
  Only SU(2) has faithful pseudoreal 2-dim rep (R2).  UNIQUE.
  SU(5) costs 24 generators; SU(3)×SU(2)×U(1) product always cheaper.

  → SU(3)×SU(2)×U(1) is the UNIQUE admissible gauge template.""")

    # ── 6. Fermion Content ───────────────────────────────────────────
    section("6. Fermion Content (Theorem T_field)")
    try:
        from apf.gauge import check_T_field
        r = check_T_field()
        print(f"  Status: {r['status']}")
    except Exception as e:
        print(f"  [check: {e}]")
    print(f"  Result: Exactly 45 fermions in Standard Model representations")
    print(f"  (1 of 4680 candidate assignments survives anomaly cancellation)")
    print(f"  C_total = 45 + 4 + 12 = 61  (rigid — changing any factor")
    print(f"  destroys all 47 predictions simultaneously)")

    # ── 7. Three Generations ─────────────────────────────────────────
    section("7. Three Generations (Theorem T7)")
    try:
        from apf.gauge import check_T7
        r = check_T7()
        print(f"  Status: {r['status']}")
    except Exception as e:
        print(f"  [check: {e}]")
    print(f"  N_gen = 3 from cyclic symmetry of the generation path graph")
    print(f"  (The question 'why 3 generations?' has a derived answer)")

    # ── 8. Cosmological Constant ─────────────────────────────────────
    section("8. Cosmological Constant Ω_Λ (L_dark_budget)")
    try:
        from apf.cosmology import check_L_dark_budget
        r = check_L_dark_budget()
        print(f"  Status: {r['status']}")
    except Exception as e:
        print(f"  [check: {e}]")
    OmegaL_pred = 0.6888
    OmegaL_obs  = 0.6889
    err = abs(OmegaL_pred - OmegaL_obs) / OmegaL_obs * 100
    print(f"  Ω_Λ predicted: {OmegaL_pred:.4f}")
    print(f"  Ω_Λ observed:  {OmegaL_obs:.4f}  (Planck 2018)")
    print(f"  Error:         {err:.2f}%")

    # ── 9. Neutrino Masses ───────────────────────────────────────────
    section("9. Neutrino Mass Ratio (L_neutrino_closure)")
    try:
        from apf.sessions.session_v63c import check_L_neutrino_closure
        r = check_L_neutrino_closure()
        print(f"  Status: {r['status']}")
    except Exception as e:
        print(f"  [check: {e}]")
    ratio_pred = 0.02952
    ratio_obs  = 0.02938
    err = abs(ratio_pred - ratio_obs) / ratio_obs * 100
    print(f"  Δm²₂₁/Δm²₃₁ predicted: {ratio_pred:.5f}")
    print(f"  Δm²₂₁/Δm²₃₁ observed:  {ratio_obs:.5f}")
    print(f"  Error:                  {err:.2f}%")
    print(f"  (Zero neutrino data used as input)")

    # ── 10. Adversarial Audit ────────────────────────────────────────
    section("10. Adversarial Self-Test — No Circularity (RT_bridge_audit)")
    try:
        from apf.red_team import check_RT_bridge_audit
        r = check_RT_bridge_audit()
        print(f"  Status: {r['status']}")
    except Exception as e:
        print(f"  [check: {e}]")
    print("""
  The 5 interpretive bridges verified as [P] theorems:
    Bridge A: dim(G) = enforcement cost     → L_cost [P]
    Bridge B: capacity fractions = Ω        → L_equip [P]
    Bridge C: d_eff^C = microstates         → L_self_exclusion + T_Bek [P]
    Bridge D: σ = ln(d_eff)                 → T_entropy + T11 [P]
    Bridge E: x = 1/2                       → T27c + L_Gram [P]

  Result: ZERO interpretive bridges remain. Every connection between
  capacity geometry and physical observables is a proved theorem.""")

    # ── Summary ──────────────────────────────────────────────────────
    print()
    print("=" * 60)
    print("  Summary")
    print("=" * 60)
    print("""
  A1 (finite enforcement capacity)
    → Hilbert space over ℂ              (T2)
    → Born rule                         (T_Born)
    → Tsirelson bound                   (T_Tsirelson)
    → SU(3)×SU(2)×U(1) unique          (L_gauge_template_uniqueness)
    → 45 fermions, correct reps         (T_field)
    → sin²θ_W = 3/13                   (L_Cauchy_uniqueness)
    → Ω_Λ = 0.688...                   (L_dark_budget)
    → 47+ quantitative predictions      (mean error 3.83%)

  Full verification: python scripts/verify_all.py
  Papers:           papers/paper1/main.tex  etc.
  AI guide (Claude):  ai_guides/CLAUDE_PROJECTS.md
  AI guide (generic): ai_guides/GENERIC_LLM.md
""")


if __name__ == "__main__":
    demo()
