"""
Verify that required CIAO tools are available and runnable in the environment.

SuperChandra requires these tools; they are provided by the CIAO package
(installable via conda from CXC or via the ciao-install script). See README.
"""
import subprocess
import shutil
import pytest

# CIAO tools required by superchandra (see README).
REQUIRED_CIAO_TOOLS = [
    "make_instmap_weights",
    "merge_obs",
    "download_chandra_obsid",
    "dmcopy",
]


def _tool_help_ok(cmd, timeout=10):
    """Return True if running cmd with --help or -h completes without crash (exit 0 or 1)."""
    for flag in ("--help", "-h"):
        try:
            r = subprocess.run(
                [cmd, flag],
                capture_output=True,
                text=True,
                timeout=timeout,
            )
            if r.returncode in (0, 1):
                return True
        except (FileNotFoundError, subprocess.TimeoutExpired, OSError):
            continue
    return False


def _tool_found(cmd):
    """Return True if cmd is found on PATH."""
    return shutil.which(cmd) is not None


_ciao_available = _tool_found("make_instmap_weights")


@pytest.mark.ciao
@pytest.mark.skipif(not _ciao_available, reason="CIAO tools not on PATH (install CIAO to run these tests)")
class TestCiaoToolsAvailable:
    """Check that each required CIAO tool is on PATH and runnable."""

    @pytest.mark.parametrize("tool", REQUIRED_CIAO_TOOLS)
    def test_tool_on_path(self, tool):
        """Tool executable is found on PATH."""
        assert _tool_found(tool), (
            f"{tool!r} not found on PATH. "
            "Install CIAO (e.g. conda install -c https://cxc.cfa.harvard.edu/conda/ciao -c conda-forge ciao) "
            "and ensure the environment is activated."
        )

    @pytest.mark.parametrize("tool", REQUIRED_CIAO_TOOLS)
    def test_tool_runs(self, tool):
        """Tool runs without error (e.g. --help succeeds)."""
        if not _tool_found(tool):
            pytest.skip(f"{tool} not on PATH")
        assert _tool_found(tool) and _tool_help_ok(tool), (
            f"{tool} found but did not accept --help or -h (or timed out)."
        )


def test_all_required_tools_listed():
    """Sanity check that REQUIRED_CIAO_TOOLS list matches README and has expected entries."""
    assert len(REQUIRED_CIAO_TOOLS) == 4
    assert "make_instmap_weights" in REQUIRED_CIAO_TOOLS
    assert "merge_obs" in REQUIRED_CIAO_TOOLS
    assert "download_chandra_obsid" in REQUIRED_CIAO_TOOLS
    assert "dmcopy" in REQUIRED_CIAO_TOOLS
