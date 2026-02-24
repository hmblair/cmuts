"""
Pytest configuration and fixtures for cmuts integration tests.
"""
from __future__ import annotations

import shutil
import sys
from pathlib import Path

import pytest

# Add tests/python to path for imports
sys.path.insert(0, str(Path(__file__).parent))


@pytest.fixture
def samtools_available() -> bool:
    """Check if samtools is available."""
    return shutil.which("samtools") is not None


@pytest.fixture
def cmuts_available() -> bool:
    """Check if cmuts binaries are available."""
    return (
        shutil.which("cmuts") is not None and
        shutil.which("_cmuts-generate-tests") is not None
    )


@pytest.fixture(autouse=True)
def require_dependencies(samtools_available: bool, cmuts_available: bool) -> None:
    """Skip tests if required dependencies are missing."""
    if not samtools_available:
        pytest.skip("samtools not found in PATH")
    if not cmuts_available:
        pytest.skip("cmuts binaries not found in PATH")
