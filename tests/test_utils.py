"""
Unit tests for standalone utility functions (validation, math, coord, extinction).
"""
import numpy as np
import pytest
from astropy.coordinates import SkyCoord

from superchandra.superchandra import (
    is_number,
    integrate_trapezoid,
    parse_coord,
    Av_to_nH,
)


class TestIsNumber:
    def test_numeric_string(self):
        assert is_number("12.5") is True
        assert is_number("-3") is True
        assert is_number("0") is True

    def test_non_numeric_string(self):
        assert is_number("abc") is False
        assert is_number("12:30:00") is False
        assert is_number("") is False

    def test_none_and_invalid(self):
        assert is_number(None) is False
        assert is_number([]) is False


class TestIntegrateTrapezoid:
    def test_simple(self):
        x = np.array([0.0, 1.0, 2.0])
        y = np.array([0.0, 1.0, 2.0])
        result = integrate_trapezoid(x, y)
        assert result == pytest.approx(2.0)  # area under y=x from 0 to 2

    def test_length_mismatch(self):
        x = np.array([0.0, 1.0])
        y = np.array([0.0, 1.0, 2.0])
        assert integrate_trapezoid(x, y) is None

    def test_two_points(self):
        x = np.array([0.0, 1.0])
        y = np.array([1.0, 1.0])
        assert integrate_trapezoid(x, y) == pytest.approx(1.0)


class TestParseCoord:
    def test_degrees(self):
        c = parse_coord(12.5, -5.2)
        assert isinstance(c, SkyCoord)
        assert c.ra.deg == pytest.approx(12.5)
        assert c.dec.deg == pytest.approx(-5.2)

    def test_hms_dms(self):
        c = parse_coord("00:50:00", "-05:12:00")
        assert isinstance(c, SkyCoord)
        assert c.ra.deg == pytest.approx(12.5)
        assert c.dec.deg == pytest.approx(-5.2)


class TestAvToNh:
    def test_conversion(self):
        # Güver & Özel: nH = 2.21e21 * Av
        assert Av_to_nH(0) == 0.0
        assert Av_to_nH(1.0) == pytest.approx(2.21e21)
        assert Av_to_nH(0.5) == pytest.approx(1.105e21)
