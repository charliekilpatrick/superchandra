# superchandra

[![Build and test](https://github.com/charliekilpatrick/superchandra/actions/workflows/build-and-test.yml/badge.svg?branch=master)](https://github.com/charliekilpatrick/superchandra/actions/workflows/build-and-test.yml)

Download and analyze Chandra/ACIS data for a given sky position: observation lookup, data download (evt2, asol, bpix, mask), merging with CIAO, and aperture photometry with configurable spectral models (e.g. blackbody, power law) for flux and upper limits.

**Author:** C. D. Kilpatrick  

---

## Citation

If you use this code in your research, please cite:

**Kilpatrick, C. D., Coulter, D. A., Dimitriadis, G., et al. 2018**, *X-ray limits on the progenitor system of the Type Ia supernova 2017ejb*, MNRAS, 481, 4123.  
- ADS: [2018MNRAS.481.4123K](https://ui.adsabs.harvard.edu/abs/2018MNRAS.481.4123K/abstract)  
- DOI: [10.1093/mnras/sty2503](https://doi.org/10.1093/mnras/sty2503)  
- arXiv: [1806.08436](https://arxiv.org/abs/1806.08436)

---

## Installation

**Requirements:** [Conda](https://docs.conda.io/en/latest/) (Miniconda or Anaconda) and Python 3.12.

1. Clone or download this repository and go to the project root (the directory containing `pyproject.toml` and `environment.yml`).

2. Create and activate the conda environment from `environment.yml`:

   ```bash
   conda env create -f environment.yml
   conda activate superchandra
   ```

3. Install the package in editable mode (dependencies are installed from `pyproject.toml`):

   ```bash
   pip install -e .
   ```

4. Verify the CLI is available:

   ```bash
   superchandra --help
   ```

The `superchandra` command runs the main pipeline. You can also run:

```bash
python -m superchandra.superchandra ra dec [options]
```

**Alternative (venv):** If you prefer a virtualenv instead of conda, use Python 3.12 to create it (`python3.12 -m venv .venv`), activate it, then run `pip install -e .` as above.

---

## Build and test (conda)

From the project root, using conda:

```bash
conda env create -f environment.yml
conda activate superchandra
pip install -e ".[test]"
pytest -m "not ciao and not network" tests/
```

To run the full test suite (including tests that require CIAO when available):

```bash
pytest -m "not network" tests/
```

To run all tests (including network-dependent ones): `pytest tests/`

---

## Dependencies (Python)

Python dependencies are declared in `pyproject.toml` and installed when you run `pip install -e .` (or `pip install -e ".[test]"` for tests):

| Package         | Purpose                          |
|----------------|----------------------------------|
| numpy          | Numerical arrays and math        |
| astropy        | Coordinates, tables, FITS, WCS   |
| requests       | HTTP requests (CDA, IRSA DUST)  |
| xmltodict      | Parse VO/XML from CDA            |
| photutils      | Aperture photometry              |
| python-dateutil| Date parsing (--before, --after) |

The `environment.yml` file defines only the conda environment name and Python version; all package dependencies are installed via `pip install -e .` from `pyproject.toml`.

---

## Requirements (external)

- **CIAO** (Chandra Interactive Analysis of Observations), version 4.11 or later, provides the tools superchandra needs:
  - `make_instmap_weights` — spectral weights for instrument maps
  - `merge_obs` — merge event lists and build merged images/expmaps
  - `download_chandra_obsid` — download evt2, asol, bpix, msk by ObsID
  - `dmcopy` — FITS manipulation (e.g. event list to image)

  You can install CIAO **via conda** (recommended) or with the **ciao-install script** from the CXC website. Both provide the same tools; conda keeps everything in one environment.

### Option A: Install CIAO with conda (recommended)

CIAO is distributed on the CXC conda channel. One environment can hold both CIAO and superchandra (install superchandra into the CIAO env with pip).

```bash
# Create env with CIAO (CIAO 4.18 uses Python 3.12; pin numpy for compatibility)
conda create -n ciao-4.18 -c https://cxc.cfa.harvard.edu/conda/ciao -c conda-forge \
  ciao sherpa ds9 ciao-contrib caldb_main marx numpy=2.3.5
conda activate ciao-4.18
pip install -e .
```

See [Installing CIAO with conda](https://cxc.cfa.harvard.edu/ciao/threads/ciao_install_conda/) for the latest version and options. After installation, the four tools above are on `PATH` when the env is activated. To verify they run:

```bash
pytest tests/test_ciao_tools.py -v
```

### Option B: Install CIAO from the CXC website

Use the [ciao-install script](http://cxc.harvard.edu/ciao/threads/ciao_install_tool/), then either put CIAO's `bin` on `PATH` and use a separate superchandra env, or install superchandra inside CIAO’s constrained environment instead:

```bash
source $CIAO_DIR/bin/ciao.bash
pip install -c $ASCDS_INSTALL/constraints.txt -e .
```

---

## Usage

Basic form: provide RA and Dec (degrees or colon-separated hours/degrees).

```bash
superchandra 12.5 -5.2
superchandra 00:50:00 -05:12:00 --download --model sss
```

---

## Command-line options

| Option | Short | Type | Default | Description |
|--------|-------|------|---------|--------------|
| `ra` | — | (positional) | — | Right ascension (degrees or HH:MM:SS). |
| `dec` | — | (positional) | — | Declination (degrees or DD:MM:SS). |
| `--makeclean` | — | flag | off | Clean up all output files from previous runs, then exit. |
| `--download` | — | flag | off | Download the raw data files for the matched observations. |
| `--before` | — | YYYY-MM-DD | — | Reject all Chandra/ACIS observations *after* this date. |
| `--after` | — | YYYY-MM-DD | — | Reject all Chandra/ACIS observations *before* this date. |
| `--clobber` | — | flag | off | Overwrite existing files when using download mode. |
| `--extinction` | — | float | 0 | Host extinction A_V (magnitudes). MW-like dust (Güver & Özel 2009) is assumed. |
| `--redshift` | `-z` | float | 0 | Host redshift for extinction. |
| `--quick` | `-q` | flag | off | Run only the first model in the chosen set (faster). |
| `--model` | `-m` | str | sss | Spectral model set: `sss` (blackbody grid), `bh` (soft/hard black hole), or `frb` (power law). |
| `--search_radius` | `-s` | float | 0.1667 | Search radius for observation lookup (degrees). |
| `--radius` | `-r` | float | (see below) | Aperture radius in arcsec for photometry (default from config, e.g. 2.216). |

For the full list and any extra options, run `superchandra --help`.

---

## Contact and support

For **bugs**, please open an [issue](https://github.com/ckilpatrick/superchandra/issues) on GitHub. For questions or feature requests, you may also email the developer:

**Charlie Kilpatrick**  
Email: [ckilpatrick@northwestern.edu](mailto:ckilpatrick@northwestern.edu)

Include your Python and CIAO versions and a short description of the issue (and, if possible, minimal steps to reproduce).

---

## License

This project is distributed under the MIT License. See the `LICENSE` file in the repository for details.
