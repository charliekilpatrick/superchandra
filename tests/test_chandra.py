"""
Unit tests for Chandra class methods.

No CIAO or network required; uses fixtures for weight file, FITS image, and expmap.
"""
import os
import argparse
import pytest

from superchandra.superchandra import Chandra


class TestChandraInit:
    def test_defaults(self, chandra):
        assert chandra.coord is None
        assert chandra.rawdir == "raw/"
        assert chandra.weightdir == "weights/"
        assert chandra.options["spectrum"]["type"] == "blackbody"
        assert len(chandra.options["spectrum"]["models"]) == 37


class TestSetSpectrum:
    def test_sss(self, chandra):
        chandra.set_spectrum("sss")
        assert chandra.options["spectrum"]["type"] == "blackbody"
        assert chandra.names[0] == "temperature"

    def test_bh(self, chandra):
        chandra.set_spectrum("bh")
        assert chandra.options["spectrum"]["type"] == "blackhole"
        assert len(chandra.options["spectrum"]["models"]) == 2

    def test_frb(self, chandra):
        chandra.set_spectrum("frb")
        assert chandra.names[0] == "gamma"
        assert len(chandra.final_phot) == 0


class TestAddOptions:
    def test_parser_creation(self, chandra):
        parser = chandra.add_options(usage="superchandra ra dec")
        assert isinstance(parser, argparse.ArgumentParser)

    def test_parser_has_expected_args(self, chandra):
        parser = chandra.add_options()
        # Check positional ra, dec and default model.
        args = parser.parse_args(["10", "20"])
        assert args.ra == "10"
        assert args.dec == "20"
        assert args.model == "sss"
        assert args.download is False


class TestMakeBanner:
    def test_prints(self, chandra, capsys):
        from superchandra.utils.logging import configure_logging
        import logging
        configure_logging(level=logging.INFO)
        chandra.make_banner("Test banner")
        out, err = capsys.readouterr()
        # Banner is logged to stderr via logging.
        combined = out + err
        assert "Test banner" in combined
        assert "#" in combined


class TestGetBolometricCorrection:
    def test_blackbody(self, chandra):
        model = {
            "temperature": 0.5,
            "weights": {"emin": 0.3, "emax": 1.0},
        }
        bccorr = chandra.get_bolometric_correction(model, spectype="blackbody")
        assert 0 < bccorr < 1

    def test_non_blackbody_returns_one(self, chandra):
        model = {"gamma": 2.0, "weights": {"emin": 0.3, "emax": 10.0}}
        assert chandra.get_bolometric_correction(model, spectype="powerlaw") == 1.0


class TestAnalyzeWeightfile:
    def test_minimal_weight_file(self, chandra, minimal_weight_file):
        avg_energy = chandra.analyze_weightfile(minimal_weight_file)
        # For linear energy 0.5--2 keV and weight 1, mean energy ~ 1.25 keV.
        assert 0.5 <= avg_energy <= 2.0


class TestDoPhotometry:
    def test_minimal_image(self, chandra, minimal_fits_image, coord_deg):
        phot, back = chandra.do_photometry(minimal_fits_image, coord_deg, radius=2.0)
        assert phot >= 0
        assert back >= 0


class TestGetExpmapValue:
    def test_minimal_expmap(self, chandra, minimal_expmap, coord_deg):
        exp = chandra.get_expmap_value(minimal_expmap, coord_deg, radius=2.0)
        assert exp > 0


class TestMergeObs:
    """merge_obs with empty file lists does not produce merged.fits."""

    def test_empty_lists_returns_false(self, chandra, temp_dir):
        workdir = temp_dir / "merge_work"
        workdir.mkdir()
        outdir = str(workdir / "merge_out") + "/"
        os.makedirs(outdir, exist_ok=True)
        weightfile = str(workdir / "weights_dummy.dat")
        with open(weightfile, "w") as f:
            f.write("0.5 1.0\n1.0 1.0\n")
        orig_cwd = os.getcwd()
        try:
            os.chdir(workdir)
            result = chandra.merge_obs(
                [], [], [], [], outdir, weightfile, coord=None, stdout=False
            )
        finally:
            os.chdir(orig_cwd)
        assert result is False
