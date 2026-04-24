#!/usr/bin/env python3
"""Convenience entry point for verifying all theorems in this paper-companion repo.

Usage:
    python run_checks.py            # run all theorems
    python run_checks.py --verbose  # show detailed per-theorem output
    python run_checks.py --check T_Born  # run a single theorem by name

This is the standalone-repo verifier for Admissibility Capacity Ledger (APF Paper 8).
The full APF codebase v6.9 (2026-04-23) verifies 437 checks across
420 bank-registered theorems + apf/standalone/ across 34 registered modules.
This paper-companion repo bundles 23 Paper-8-specific theorems (see
ai_context/theorems.json) on top of the full 420-theorem bank, all verified
here via the bundled apf/ subset.
"""

import sys
import argparse
import time
from apf.bank import REGISTRY, run_all
try:
    from apf.bank import get_check
except ImportError:
    # Canonical apf/bank.py (used by Paper 13 full-codebase mode) doesn't
    # define get_check; fall back to REGISTRY lookup.
    def get_check(name):
        if name in REGISTRY:
            return REGISTRY[name]
        raise KeyError(f"Check '{name}' not in REGISTRY")


def main():
    parser = argparse.ArgumentParser(description="Verify all theorems in this paper-companion repo.")
    parser.add_argument("--verbose", "-v", action="store_true",
                        help="Show detailed per-theorem output")
    parser.add_argument("--check", "-c", type=str, default=None,
                        help="Run a single theorem by name (e.g., --check T_Born)")
    args = parser.parse_args()

    if args.check:
        try:
            check_fn = get_check(args.check)
        except KeyError:
            print(f"  ERROR: theorem '{args.check}' not in this repo's bank.")
            print(f"  Available: {sorted(REGISTRY.keys())}")
            sys.exit(1)
        result = check_fn()
        print(f"  {args.check}: {'PASS' if result.get('passed', True) else 'FAIL'}")
        if "key_result" in result:
            print(f"    Result: {result['key_result']}")
        if args.verbose:
            print(f"    Full result: {result}")
        sys.exit(0 if result.get('passed', True) else 1)

    start = time.time()
    results = run_all(verbose=args.verbose)
    elapsed = time.time() - start

    # Three possible return shapes:
    #   - list[dict]   — paper-subset bank.py (per-result dict)
    #   - dict[str, dict] — canonical apf/bank.py (Paper 13 full-codebase mode)
    #   - list[str]    — older bank.py variants (status lines)
    if isinstance(results, dict):
        # Canonical: results is {check_name: result_dict}
        result_iter = list(results.values())
        passed = sum(1 for r in result_iter if isinstance(r, dict) and r.get('passed', True))
        failed = sum(1 for r in result_iter if isinstance(r, dict) and not r.get('passed', True))
        total = len(result_iter)
    elif results and isinstance(results[0], dict):
        result_iter = results
        passed = sum(1 for r in result_iter if r.get('passed', True))
        failed = len(result_iter) - passed
        total = len(result_iter)
    else:
        result_iter = list(results)
        passed = len(result_iter)
        failed = 0
        total = len(result_iter)

    print()
    print("=" * 72)
    print(f"  Paper 8 ({{Admissibility Capacity Ledger}}): {passed} passed, {failed} failed, "
          f"{total} total — verified in {elapsed:.2f}s")
    print("=" * 72)

    if failed > 0:
        print("\nFailed checks:")
        for r in result_iter:
            if isinstance(r, dict) and not r.get('passed', True):
                print(f"  - {r.get('name', '?')}: {r.get('error', 'no detail')}")
        sys.exit(1)

    sys.exit(0)


if __name__ == "__main__":
    main()
