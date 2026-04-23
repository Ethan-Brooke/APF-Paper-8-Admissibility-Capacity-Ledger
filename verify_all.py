#!/usr/bin/env python3
"""APF Full Verification Script.

Runs all check_ functions across all registered modules and prints a
complete scorecard with prediction accuracy summary.

Usage:
    python scripts/verify_all.py
    python scripts/verify_all.py --module core
    python scripts/verify_all.py --verbose
    python scripts/verify_all.py --no-scorecard
"""

import sys
import time
import argparse
import importlib
import traceback
from collections import defaultdict

# ── Module registry ──────────────────────────────────────────────────────────
# Mirrors the load order in apf/bank.py, plus standalone files.

MODULES = [
    # Core framework (Paper 1)
    ("apf.core",                       "core.py"),
    # Standard Model sector (Paper 2)
    ("apf.gauge",                      "gauge.py"),
    ("apf.generations",                "generations.py"),
    ("apf.spacetime",                  "spacetime.py"),
    ("apf.gravity",                    "gravity.py"),
    ("apf.plec",                       "plec.py"),
    ("apf.unification",                "unification.py"),
    ("apf.unification_three_levels",   "unification_three_levels.py"),
    ("apf.unification_projection_essentiality", "unification_projection_essentiality.py"),
    ("apf.subspace_functors",          "subspace_functors.py"),
    ("apf.fractional_reading",         "fractional_reading.py"),
    ("apf.lambda_absolute",            "lambda_absolute.py"),
    ("apf.lambda_operator_derivation", "lambda_operator_derivation.py"),
    ("apf.thermal_absolute",           "thermal_absolute.py"),
    ("apf.phase_14d3_completions",     "phase_14d3_completions.py"),
    ("apf.quantum_operator_derivation","quantum_operator_derivation.py"),
    ("apf.i4_composition",             "i4_composition.py"),
    ("apf.horizon_joint_bridge",       "horizon_joint_bridge.py"),
    ("apf.crystal",                    "crystal.py"),
    ("apf.crystal_metrics",            "crystal_metrics.py"),
    ("apf.killed_rivals",              "killed_rivals.py"),
    ("apf.cosmology",                  "cosmology.py"),
    # Extended checks
    ("apf.supplements",                "supplements.py"),
    ("apf.majorana",                   "majorana.py"),
    ("apf.internalization",            "internalization.py"),
    ("apf.internalization_geo",        "internalization_geo.py"),
    ("apf.extensions",                 "extensions.py"),
    ("apf.validation",                 "validation.py"),
    # Adversarial
    ("apf.red_team",                   "red_team.py"),
    # Sessions (incremental development)
    ("apf.session_v63c",               "session_v63c.py"),
    ("apf.session_qg",                 "session_qg.py"),
    ("apf.session_nnlo",               "session_nnlo.py"),
    ("apf.session_delta_pmns",         "session_delta_pmns.py"),
    ("apf.session_cosmo_update",       "session_cosmo_update.py"),
    ("apf.session_phase2_confrontation","session_phase2_confrontation.py"),
    # Standalone lemmas
    ("apf.standalone.L_Cauchy_uniqueness",    "standalone/L_Cauchy_uniqueness.py"),
    ("apf.standalone.L_CKM_resolution_limit", "standalone/L_CKM_resolution_limit.py"),
    ("apf.standalone.phase1_seesaw_closure",  "standalone/phase1_seesaw_closure.py"),
    ("apf.standalone.phase5_theorem_R_audit", "standalone/phase5_theorem_R_audit.py"),
]

# Modules whose checks populate the DAG cache. When --module filter is used,
# these run as a silent prelude (if not in the filtered set) so that
# downstream DAG consumers see a populated DAG. Derived from the 2026-04-20
# DAG-key audit: these are exactly the modules containing dag_put calls.
PRELUDE_MODULES = {
    "apf.core",
    "apf.gauge",
    "apf.spacetime",
    "apf.gravity",
    "apf.generations",
    "apf.cosmology",
    "apf.extensions",
    "apf.supplements",
}

# ── Helpers ──────────────────────────────────────────────────────────────────

def classify(name):
    n = name.replace("check_", "")
    if n.startswith("T_") or (n.startswith("T") and n[1:2].isdigit()):
        return "theorem"
    if n.startswith("L_"):
        return "lemma"
    if n.startswith("RT_"):
        return "red_team"
    return "other"


def run_module(mod_name, verbose=False):
    results = []
    try:
        mod = importlib.import_module(mod_name)
    except ModuleNotFoundError as e:
        return results, f"SKIP ({e})"
    except Exception as e:
        return results, f"IMPORT ERROR: {e}"

    checks = sorted(
        (name, fn)
        for name, fn in vars(mod).items()
        if name.startswith("check_") and callable(fn)
    )

    for name, fn in checks:
        status = "PASS"
        error_msg = None
        try:
            fn()
        except Exception as e:
            status = "FAIL"
            error_msg = str(e)
            if verbose:
                traceback.print_exc()
        results.append({
            "name": name,
            "category": classify(name),
            "status": status,
            "error": error_msg,
        })

    return results, None


def run_prelude_silently(prelude_mod_names):
    """Run DAG-producing modules silently to populate the DAG cache.

    Used when `--module` filter is active and the filtered set excludes
    upstream DAG producers. Runs each prelude module's check_* functions
    best-effort: exceptions inside checks don't propagate (the DAG will
    simply be missing the keys that the failing check would have
    written, which will produce a clear error downstream via the
    dag_has guards in consumer code).

    Returns
    -------
    tuple
        (n_modules_run, n_checks_run) for a one-line banner.
    """
    n_modules = 0
    n_checks = 0
    for mod_name in prelude_mod_names:
        try:
            mod = importlib.import_module(mod_name)
        except Exception:
            continue
        n_modules += 1
        check_names = sorted(n for n in dir(mod) if n.startswith('check_'))
        for name in check_names:
            try:
                fn = getattr(mod, name)
                if callable(fn):
                    fn()
                    n_checks += 1
            except Exception:
                pass  # best-effort — missing DAG keys will surface in the filtered run
    return n_modules, n_checks


def print_prediction_scorecard():
    try:
        from apf.validation import check_L_prediction_catalog
        r = check_L_prediction_catalog()
        print()
        print("Prediction scorecard")
        print("-" * 44)
        if isinstance(r, dict) and "summary" in r:
            print(r["summary"])
        else:
            print("  (run apf/validation.py directly for full catalog)")
    except Exception as e:
        print(f"\nPrediction scorecard unavailable: {e}")


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="APF Full Verification")
    parser.add_argument("--module", help="Run only checks matching this string (e.g. 'core')")
    parser.add_argument("--verbose", "-v", action="store_true", help="Print full tracebacks on failure")
    parser.add_argument("--no-scorecard", action="store_true", help="Skip prediction scorecard")
    args = parser.parse_args()

    print("APF v6.9 — Full Verification")
    print("=" * 55)

    modules_to_run = MODULES
    if args.module:
        modules_to_run = [(m, f) for m, f in MODULES if args.module.lower() in m.lower()]
        if not modules_to_run:
            print(f"No modules matching '{args.module}' found.")
            sys.exit(1)

        # If the filter excludes upstream DAG-producing modules, silently
        # pre-run them so downstream consumers see a populated DAG. Without
        # this, e.g. `--module unification` fails because L_count in
        # apf.gauge has not yet written C_total, and the 2026-04-20 dag_has
        # guards in acc_SM raise loudly on missing keys (by design).
        filtered_names = {m for m, _ in modules_to_run}
        missing_prelude = [m for m in MODULES
                           if m[0] in PRELUDE_MODULES
                           and m[0] not in filtered_names]
        if missing_prelude:
            prelude_names = [m for m, _ in missing_prelude]
            n_mods, n_checks = run_prelude_silently(prelude_names)
            print(f"[prelude: ran {n_checks} checks across {n_mods} "
                  f"upstream DAG producer(s) to populate cache for filtered run]")
            print()

    all_results = []
    module_summaries = []
    t_start = time.time()

    for mod_name, file_label in modules_to_run:
        results, err = run_module(mod_name, verbose=args.verbose)
        n_pass = sum(1 for r in results if r["status"] == "PASS")
        n_fail = sum(1 for r in results if r["status"] == "FAIL")
        n = len(results)

        if err:
            status_str = err
        elif n_fail == 0 and n > 0:
            status_str = "PASS"
        elif n == 0:
            status_str = "no checks"
        else:
            status_str = f"FAIL ({n_fail}/{n})"

        module_summaries.append((file_label, n, status_str))
        all_results.extend(results)

        if n_fail > 0:
            for r in results:
                if r["status"] == "FAIL":
                    print(f"  FAIL  {r['name']}")
                    if r["error"]:
                        print(f"        {r['error'][:120]}")

    elapsed = time.time() - t_start

    print(f"\n{'Module':<47} {'Checks':>6}  {'Status'}")
    print(f"{'──────':<47} {'──────':>6}  {'──────'}")
    for label, n, status in module_summaries:
        print(f"{label:<47} {n:>6}  {status}")

    cats = defaultdict(lambda: {"pass": 0, "fail": 0})
    for r in all_results:
        cats[r["category"]]["pass" if r["status"] == "PASS" else "fail"] += 1

    total_pass = sum(r["status"] == "PASS" for r in all_results)
    total_fail = sum(r["status"] == "FAIL" for r in all_results)
    total = len(all_results)

    print()
    print("=" * 55)
    print(f"Total checks:    {total:>4}")
    print(f"  Theorems  [T]: {cats['theorem']['pass']:>4}  pass   {cats['theorem']['fail']:>3}  fail")
    print(f"  Lemmas    [L]: {cats['lemma']['pass']:>4}  pass   {cats['lemma']['fail']:>3}  fail")
    print(f"  Red-team [RT]: {cats['red_team']['pass']:>4}  pass   {cats['red_team']['fail']:>3}  fail")
    print(f"  Other:         {cats['other']['pass']:>4}  pass   {cats['other']['fail']:>3}  fail")
    print()

    if total_fail == 0 and total > 0:
        print(f"✓ All {total} checks PASSED  ({elapsed:.1f}s)")
    elif total == 0:
        print("No checks found — check module paths.")
    else:
        print(f"✗ FAILURES: {total_fail}/{total} checks FAILED  ({elapsed:.1f}s)")
        print("\nFailed checks:")
        for r in all_results:
            if r["status"] == "FAIL":
                print(f"  {r['name']}")
                if r["error"]:
                    print(f"    {r['error'][:200]}")

    if not args.no_scorecard:
        print_prediction_scorecard()

    sys.exit(0 if total_fail == 0 else 1)


if __name__ == "__main__":
    main()
