"""
Tests for the superchandra CLI (main entry point: --help, required ra/dec, importable main).
"""
import subprocess
import sys
import pytest


def test_superchandra_help():
    """CLI --help exits 0 and prints usage."""
    result = subprocess.run(
        [sys.executable, "-m", "superchandra.superchandra", "--help"],
        capture_output=True,
        text=True,
        timeout=10,
    )
    assert result.returncode == 0
    assert "ra" in result.stdout or "RA" in result.stdout
    assert "dec" in result.stdout or "DEC" in result.stdout


def test_superchandra_requires_ra_dec():
    """CLI without ra/dec exits non-zero and prints usage."""
    result = subprocess.run(
        [sys.executable, "-m", "superchandra.superchandra"],
        capture_output=True,
        text=True,
        timeout=10,
    )
    assert result.returncode != 0
    assert "superchandra" in (result.stdout + result.stderr).lower() or "ra" in (result.stdout + result.stderr).lower()


def test_main_importable():
    """main() is callable from the package."""
    from superchandra.superchandra import main
    assert callable(main)
