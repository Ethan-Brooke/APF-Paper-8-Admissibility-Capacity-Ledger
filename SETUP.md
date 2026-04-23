# Setup and Verification

## Requirements

- Python 3.8 or later
- **numpy ≥ 1.20** — linear algebra, array operations; used across most of the bank
- **scipy ≥ 1.7** — special functions, two-loop RG integrals, spectral-action moments; specifically required by 4 checks (`L_SA_moments`, `L_ST_index`, `L_mc_mt_twoloop_RG`, `L_spectral_action_coefficients`)
- Python stdlib (`fractions`, `math`, `itertools`, `functools`) for the exact-rational core

`pip install -e .` reads dependencies from `setup.py` automatically. As of v6.8 (frozen 2026-04-18), the engine is **not** stdlib-only — earlier setup notes describing it as such were stale. The numerical components were added incrementally; the canonicalization pass surfaced this and updated the contract.

## Installation

```bash
# Clone the repository
git clone https://github.com/Ethan-Brooke/APF-Complete.git
cd APF-Complete

# Install the apf package in editable mode
pip install -e .

# Verify the install
python -c "import apf; print('APF package loaded successfully')"
```

## Package Structure

The `apf/` directory is a standard Python package. After `pip install -e .` it is importable from anywhere as `import apf` or `from apf.core import check_T2`.

A minimal `setup.py` / `pyproject.toml` is included. If you prefer not to install, you can also run scripts directly from the repo root by ensuring the root is on your `PYTHONPATH`:

```bash
export PYTHONPATH=/path/to/APF-Complete:$PYTHONPATH
python scripts/verify_all.py
```

## Running All Checks

```bash
python scripts/verify_all.py
```

This runs all 349 check functions across all modules and prints a full scorecard. Expected output:

```
APF v6.7 — Full Verification
===========================================
Module                Checks   Status
------                ------   ------
core.py                   48   PASS
gauge.py                  29   PASS
cosmology.py              17   PASS
gravity.py                 9   PASS
spacetime.py               8   PASS
majorana.py               10   PASS
generations__8_.py        86   PASS
extensions.py              7   PASS
supplements.py            73   PASS
validation.py              9   PASS
red_team.py               19   PASS
internalization.py         3   PASS
internalization_geo.py     4   PASS
sessions (6 files)        27   PASS
===========================================
Total checks:            349   ALL PASS

Prediction scorecard
--------------------
Tested predictions:   39
Within 1%:            24   (61.5%)
Within 5%:            36   (92.3%)
Mean absolute error:  3.83%
Median error:         0.37%
```

## Running Individual Modules

```python
# Run the Paper 1 core checks
from apf.core import *
result = check_T2()
print(result['name'], '—', result['status'])

# Run the gauge sector
from apf.gauge import check_L_gauge_template_uniqueness
r = check_L_gauge_template_uniqueness()
print(r['summary'])

# Run all red-team adversarial tests
from apf.red_team import *
check_RT_bridge_audit()
check_RT_adversarial_Ngen4()
check_RT_R3_no_U1()
```

## Running the 5-Minute Demo

```bash
python scripts/quick_demo.py
```

This prints the 10 headline results with their derivation chains and experimental comparisons — a good first orientation to what APF actually produces.

## LaTeX Papers

The `papers/` directory contains `.tex` source files. To compile:

```bash
cd papers/paper1
pdflatex -interaction=nonstopmode main.tex
pdflatex -interaction=nonstopmode main.tex
pdflatex -interaction=nonstopmode main.tex
```

Three passes are needed for cross-references and page numbers to stabilize. The supplement is a standalone document compiled the same way.

## Troubleshooting

**Import error: `No module named 'apf'`**  
Run `pip install -e .` from the repo root, or add the repo root to `PYTHONPATH`.

**`CheckFailure` raised during verify_all**  
A theorem check has failed. The error message names the specific check and the assertion that failed. This should not happen on a clean clone — if it does, please report it.

**Slow runtime**  
The full suite runs in under 60 seconds on any modern laptop. The `generations__8_.py` module (86 checks, Gram matrix computations) takes the longest.
