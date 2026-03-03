"""
Integration tests that require CIAO (e.g. make_instmap_weights).

Skipped when CIAO is not installed or not on PATH.
"""
import os
import shutil
import pytest

from superchandra.superchandra import Chandra


def ciao_available():
    """Return True if make_instmap_weights is on PATH (CIAO environment is set up)."""
    return shutil.which("make_instmap_weights") is not None


requires_ciao = pytest.mark.skipif(
    not ciao_available(),
    reason="CIAO tools (make_instmap_weights) not on PATH",
)


@requires_ciao
@pytest.mark.ciao
class TestMakeInstmapWeightsWithCiao:
    """Test make_instmap_weights when CIAO is available."""

    def test_blackbody_weight_file_created(self, chandra, temp_dir):
        weightfile = str(temp_dir / "weights_t0.5.nh0.02.dat")
        model = {
            "name": "t0.5000",
            "temperature": 0.5,
            "weights": {"emin": 0.3, "emax": 1.0, "ewidth": 0.02},
        }
        galnh = 0.02
        hostnh = 0.0
        z = 0.0
        result = chandra.make_instmap_weights(
            weightfile, model, galnh, hostnh, z, stdout=False
        )
        assert result is True
        assert os.path.isfile(weightfile)
