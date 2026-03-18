"""
Pytest configuration and fixtures for cmuts integration tests.
"""

from __future__ import annotations

import shutil
import sys
from pathlib import Path

import pytest

# Add tests/python and src/python to path for imports
sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src" / "python"))


@pytest.fixture
def samtools_available() -> bool:
    """Check if samtools is available."""
    return shutil.which("samtools") is not None


@pytest.fixture
def cmuts_available() -> bool:
    """Check if cmuts binaries are available."""
    return shutil.which("cmuts") is not None and shutil.which("_cmuts-generate-tests") is not None


@pytest.fixture(autouse=True)
def require_dependencies(
    request: pytest.FixtureRequest,
    samtools_available: bool,
    cmuts_available: bool,
) -> None:
    """Skip external-tool-dependent tests if required binaries are missing."""
    if request.node.get_closest_marker("no_external_dependencies") is not None:
        return

    if not samtools_available:
        pytest.skip("samtools not found in PATH")
    if not cmuts_available:
        pytest.skip("cmuts binaries not found in PATH")
