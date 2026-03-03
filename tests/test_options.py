"""
Tests for superchandra.utils.options (CLI option definitions and get_parser behavior).
"""
import argparse
import pytest

from superchandra.utils.options import get_parser


class TestGetParser:
    """get_parser() returns a parser with all superchandra options."""

    def test_returns_parser(self):
        p = get_parser()
        assert isinstance(p, argparse.ArgumentParser)

    def test_accepts_usage(self):
        p = get_parser(usage="superchandra ra dec")
        assert p is not None

    def test_extends_existing_parser(self):
        base = argparse.ArgumentParser()
        base.add_argument("--extra", action="store_true")
        p = get_parser(parser=base)
        args = p.parse_args(["1.0", "2.0", "--extra"])
        assert args.ra == "1.0"
        assert args.dec == "2.0"
        assert args.extra is True

    def test_positional_ra_dec(self):
        p = get_parser()
        args = p.parse_args(["12.5", "-5.2"])
        assert args.ra == "12.5"
        assert args.dec == "-5.2"

    def test_defaults(self):
        p = get_parser()
        args = p.parse_args(["0", "0"])
        assert args.download is False
        assert args.clobber is False
        assert args.model == "sss"
        assert args.quick is False
        assert args.extinction == 0.0
        assert args.redshift == 0.0
        assert args.search_radius == 0.1667
        assert args.radius is None
        assert args.before is None
        assert args.after is None
        assert args.makeclean is False

    def test_short_flags(self):
        p = get_parser()
        args = p.parse_args(["0", "0", "-m", "bh", "-z", "0.5", "-q", "-s", "0.2", "-r", "3.0"])
        assert args.model == "bh"
        assert args.redshift == 0.5
        assert args.quick is True
        assert args.search_radius == 0.2
        assert args.radius == 3.0

    def test_flags_store_true(self):
        p = get_parser()
        args = p.parse_args(["0", "0", "--download", "--clobber", "--makeclean"])
        assert args.download is True
        assert args.clobber is True
        assert args.makeclean is True

    def test_before_after_dates(self):
        p = get_parser()
        args = p.parse_args(["0", "0", "--before", "2020-01-01", "--after", "2019-06-01"])
        assert args.before == "2020-01-01"
        assert args.after == "2019-06-01"

    def test_extinction(self):
        p = get_parser()
        args = p.parse_args(["0", "0", "--extinction", "0.5"])
        assert args.extinction == 0.5

    def test_model_values(self):
        p = get_parser()
        for model in ("sss", "bh", "frb"):
            args = p.parse_args(["0", "0", "--model", model])
            assert args.model == model
