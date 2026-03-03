# Dependency validation

This document records how `pyproject.toml` dependencies align with the installable package, unit tests, and recommended environment.

## pyproject.toml dependencies

| Dependency         | Version (pyproject) | Used by |
|--------------------|---------------------|---------|
| numpy              | >=1.22              | superchandra.py, utils/math.py, utils/spectrum.py |
| astropy            | >=5.3               | superchandra.py (votable, coordinates, table, time, fits, wcs), utils/coord.py, tests (conftest, test_utils) |
| requests           | >=2.28              | superchandra.py (get_obstable), utils/extinction.py (get_MW_Av) |
| xmltodict          | >=0.12              | superchandra.py (VO table parse), utils/extinction.py (IRSA DUST) |
| photutils          | >=1.0               | superchandra.py (SkyCircularAperture, SkyCircularAnnulus, aperture_photometry) |
| python-dateutil    | >=2.8                | superchandra.py (dateparse for --before, --after) |
| pytest             | >=7.0 (optional [test]) | All tests |

Stdlib-only imports (not in pyproject): `os`, `sys`, `shutil`, `glob`, `tempfile`, `warnings`, `argparse`, `subprocess`, `tempfile`.

## Consistency checks

- **Package imports:** Every third-party import in `superchandra` and `superchandra.utils` is declared in `pyproject.toml` (core or test).
- **Unit tests:** Tests use only pytest (optional), main dependencies (numpy, astropy), and stdlib. No missing test dependencies.
- **Python version:** `requires-python = ">=3.12,<3.13"` matches `environment.yml` (python=3.12) and README (Python 3.12). It also aligns with CIAO 4.18 (Python 3.12).

## Recommended environment

- **Conda (no CIAO):** `environment.yml` → Python 3.12, then `pip install -e ".[test]"`.
- **Conda (with CIAO 4.18):** See README Option A; use the CIAO env’s Python 3.12 and pin `numpy=2.3.5` there (numpy 2.4+ incompatible with CIAO).

## Version bounds (current)

Lower bounds are conservative and compatible with current photutils and astropy. No upper bounds are set in pyproject; when installing into a CIAO conda env, pin numpy in that env as documented in the README.
