# Test suite

Recommended: use the conda environment, then install with test extras (see main README “Build and test”):

```bash
conda activate superchandra
pip install -e ".[test]"
```

- **Unit tests** (no external services): `pytest -m "not ciao and not network" tests/`
- **Verify CIAO tools** (after installing CIAO via conda or ciao-install): `pytest tests/test_ciao_tools.py -v`
- **All tests except network**: `pytest -m "not network" tests/` (requires CIAO on PATH for `ciao`-marked tests)
- **Full suite** (including network): `pytest tests/` (requires network for CDA/IRSA)

Markers:
- `ciao`: needs CIAO tools (e.g. `make_instmap_weights`)
- `network`: needs network (CDA, IRSA)
