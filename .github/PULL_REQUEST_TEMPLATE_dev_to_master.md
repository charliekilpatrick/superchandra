# Pull Request: dev → master

## Summary

Updates from the `dev` branch: improve repo hygiene (ignore virtualenv dirs), fix build and test warnings (license declaration, numerical stability in bolometric correction).

## Changes

### .gitignore
- **Ignore `.venv*` directories** (e.g. `.venv312`) so local virtualenvs are not tracked. The pattern `.venv*/` was added alongside the existing `.venv/` entry.

### pyproject.toml
- **License declaration:** Replaced deprecated TOML table form `license = { text = "MIT" }` with the SPDX string `license = "MIT"`. This removes the SetuptoolsDeprecationWarning and aligns with the recommended format (setuptools ≥77 and packaging guidelines).

### superchandra/utils/spectrum.py
- **Bolometric correction:** In `get_bolometric_correction`, the Planck integrand can overflow in `np.exp()` for small wavelength×temperature. The computation is now wrapped in `np.errstate(over='ignore', invalid='ignore')` so overflow/invalid in the exponential is suppressed. The integrated result is unchanged (overflow only occurs in negligible tails). This removes the RuntimeWarning during tests (e.g. `TestGetBolometricCorrection::test_blackbody`).

## Testing

- [x] Build: `python -m build --outdir dist` succeeds with no license deprecation warning.
- [x] Unit tests: `pytest -m "not ciao and not network" tests/` — 46 passed, no warnings.
- [x] No new linter issues.

## Checklist

- [x] Changes are limited to dev-only fixes (gitignore, build warning, test warning).
- [x] No breaking API or behavior changes.
- [x] Python 3.12 required; consistent with existing `requires-python` and README.
