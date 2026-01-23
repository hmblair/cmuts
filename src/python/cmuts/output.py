"""Output formatting and statistics display.

This module provides formatting utilities for CLI output, including
ANSI terminal formatting and statistics display for probing data.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, SupportsFloat, Union

import numpy as np

if TYPE_CHECKING:
    from .internal import ProbingData

__all__ = [
    "BOLD",
    "RESET",
    "div",
    "space",
    "format_s",
    "format_i",
    "format_f",
    "title",
    "subtitle",
    "stats",
]


# ANSI escape codes for terminal formatting

BOLD = "\033[1m"
RESET = "\033[0m"

# Formatting constants

SPACE = " "
DIV = "─"

# Heatmap index for deletion rate
IX_DEL = 4


def div() -> str:
    """Return a divider line for CLI output."""
    return SPACE * 6 + DIV * 35


def space(n: int) -> str:
    """Return n spaces."""
    return SPACE * n


def format_s(seq: str, offset: int = 8) -> str:
    """Format a string with offset spacing."""
    return space(offset) + seq + "\n"


def format_i(seq: str, value: SupportsFloat, offset: int = 8, width: int = 40) -> str:
    """Format an integer value with label and alignment."""
    v = float(value)
    gap = width - len(f"{int(v):,}") - len(seq + ":") - offset
    return space(offset) + seq + ":" + space(gap) + f"{int(v):,}" + "\n"


def format_f(
    seq: str, value: SupportsFloat, offset: int = 8, width: int = 40, prec: int = 2
) -> str:
    """Format a float value with label, alignment, and precision."""
    v = float(value)
    gap = width - len(f"{v:.{prec}f}") - len(seq + ":") - offset
    return space(offset) + seq + ":" + space(gap) + f"{v:.{prec}f}" + "\n"


def title(name: str, version: str) -> str:
    """Format a title banner with name and version."""
    return space(8) + f"{name} version {version}\n" + div()


def subtitle(name: str) -> str:
    """Format a subtitle for statistics output."""
    if name:
        return space(8) + f"{BOLD}Statistics for {name}:{RESET}"
    else:
        return space(8) + f"{BOLD}Statistics:{RESET}"


def _print_stats_single(
    data: "ProbingData",
) -> str:
    """Format statistics for a single ProbingData object."""
    multi = data.reactivity.shape[0] > 1

    # Reads
    mr = np.mean(data.reads)
    medr = np.median(data.reads)
    dropout = (data.reads == 0).mean()

    # SNR
    CUTOFF = 1
    msnr = np.mean(data.snr)
    fsnr = np.mean(data.snr > CUTOFF)

    # Reactivity
    mrr = np.mean(data.reactivity[data.mask])

    # Error
    me = np.mean(data.error[data.mask])

    # Deletions
    del_r = np.mean(data.heatmap[:, IX_DEL])

    # Mean-to-median
    with np.errstate(divide="ignore", invalid="ignore"):
        med2m = mr / medr

    out = ""

    if multi:
        out += format_f("Mean reads", mr, 10)
        out += format_i("Median reads", medr, 10)
        out += format_f("Mean-to-Median", med2m, 10)
        out += format_f("Mean reactivity", mrr, 10, prec=3)
        out += format_f("Mean error", me, 10, prec=3)
        out += format_f("Deletion rate", del_r, 10, prec=3)
        out += format_f("Dropout fraction", dropout, 10)
        out += format_f("Mean SNR", msnr, 10)
        out += format_f("SNR > 1", fsnr, 10)
        if data.pairwise_snr is not None:
            mpsnr = np.mean(data.pairwise_snr)
            fpsnr = np.mean(data.pairwise_snr > CUTOFF)
            out += format_f("Mean pairwise SNR", mpsnr, 10)
            out += format_f("Pairwise SNR > 1", fpsnr, 10)
    else:
        out += format_i("Reads", mr, 10)
        out += format_f("Mean reactivity", mrr, 10, prec=3)
        out += format_f("Mean error", me, 10, prec=3)
        out += format_f("Deletion rate", del_r, 10, prec=3)
        out += format_f("SNR", msnr, 10)
        if data.pairwise_snr is not None:
            out += format_f("Pairwise SNR", np.mean(data.pairwise_snr), 10)

    return out


def stats(
    mod: "ProbingData",
    nomod: Union["ProbingData", None],
    combined: "ProbingData",
) -> None:
    """Print statistics for modified, unmodified, and combined probing data."""
    multi = combined.reactivity.shape[0] > 1

    nseqs = combined.reads.shape[0]
    seqlen = combined.reactivity.shape[1]
    reads = combined.reads.sum()

    print(format_i("References", nseqs, 8), end="")
    print(format_i("Reference length", seqlen, 8), end="")
    if multi:
        print(format_i("Total reads", reads, 8), end="")

    print(format_s(f"{BOLD}Treated (MOD):{RESET}", 8), end="")
    print(_print_stats_single(mod), end="")

    if nomod is not None:
        print(format_s(f"{BOLD}Untreated (NOMOD):{RESET}", 8), end="")
        print(_print_stats_single(nomod), end="")
        print(format_s(f"{BOLD}Combined:{RESET}", 8), end="")
        print(_print_stats_single(combined), end="")

    print(div())
    print()
