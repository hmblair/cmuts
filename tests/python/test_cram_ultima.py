"""
Regression test for decoding real Ultima Genomics CRAM files.

`tests/python/fixtures/ultima_softclip.cram` is a 6 KB *lossless* slice of a
real Ultima `sorter`-encoded CRAM: the first data container copied verbatim
(original codecs preserved — HUFFMAN/SUBEXP/BETA, gzip/rans4x8) with the header
trimmed to a single reference. The reads are long, mostly soft-clipped (only a
~40-60 nt core aligns to the 207 nt reference).

cmuts's hand-rolled CRAM decoder mis-decodes these records (the FN/RL data
series are absent for these records and the documented fallback is not
implemented), yielding a corrupt CIGAR whose reference span runs past the
reference, so `cmuts core` aborts the entire file. The synthetic
`cmuts-generate` CRAMs only ever emit MATCH/MISMATCH/INS/DEL with all series
present, so they never exercise this path.

This test currently FAILS and should pass once the decoder is fixed.
"""
from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import h5py
import numpy as np
import pytest

# Needs only the cmuts CLI, not samtools (the fixture is pre-built).
pytestmark = pytest.mark.no_external_dependencies

_CMUTS = shutil.which("cmuts")
FIXTURES = Path(__file__).parent / "fixtures"


@pytest.fixture(autouse=True)
def _require_cmuts() -> None:
    if _CMUTS is None:
        pytest.skip("cmuts not found in PATH")


def test_ultima_softclip_cram_decodes(tmp_path: Path) -> None:
    # Copy the fixture so the generated .cmfa index lands in tmp, not the repo.
    cram = tmp_path / "aln.cram"
    fasta = tmp_path / "ref.fasta"
    shutil.copy(FIXTURES / "ultima_softclip.cram", cram)
    shutil.copy(FIXTURES / "ultima_softclip.fasta", fasta)
    out = tmp_path / "out.h5"

    proc = subprocess.run(
        ["cmuts", "core", "-f", str(fasta), "-o", str(out), "--overwrite",
         "--min-mapq", "0", str(cram)],
        capture_output=True, text=True,
    )
    log = proc.stdout + proc.stderr

    # The decoder must process the file rather than aborting it on the first
    # mis-decoded record.
    assert "Error processing the file" not in log, log
    assert "exceeds reference sequence length" not in log, log
    assert "Invalid read length" not in log, log
    assert "0 of 1 files were processed" not in log, log

    # And the aligned cores must actually be counted.
    with h5py.File(out, "r") as f:
        found: list[str] = []
        f.visititems(lambda n, o: found.append(n) if n.endswith("counts-1d") else None)
        assert found, "no counts-1d dataset written"
        total = float(np.asarray(f[found[0]][:]).sum())
    assert total > 0, "file decoded to zero counts (records were dropped)"
