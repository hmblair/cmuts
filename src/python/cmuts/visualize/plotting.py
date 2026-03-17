import os
from typing import Any, Union

import dask.array as da
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm, SymLogNorm

from cmuts.internal import ProbingData

# Type alias for arrays that can be either numpy or dask
ArrayType = Union[np.ndarray, da.Array]

# Figures directory


FIGURES = "figures"


# Font and color for plotting


plt.rcParams["font.family"] = "Helvetica"
LABEL_SIZE = 14
TITLE_SIZE = 15
TICK_SIZE = 13
LEGEND_SIZE = 12
DPI = 400

CMAP = "RdPu"
CMAPS = [
    plt.get_cmap("RdPu"),
    plt.get_cmap("PiYG"),
]


def _line_plot(data: ArrayType, cmap: Any = None, label: str = "") -> None:
    if cmap is None:
        cmap = plt.get_cmap("RdPu")
    data = np.asarray(data)
    fill = cmap(0.3)
    draw = cmap(0.8)

    n = data.shape[0]
    x = np.arange(n)

    plt.fill_between(x, data, alpha=0.5, color=fill)
    if label:
        plt.plot(data, color=draw, linewidth=1, label=label)
    else:
        plt.plot(data, color=draw, linewidth=1)

    plt.grid(axis="y", alpha=0.5)
    plt.tick_params(axis="both", labelsize=TICK_SIZE)
    plt.xlim(0, n)


def _matrix_plot(data: np.ndarray, norm, label: str, cmap: str = CMAP) -> None:
    im = plt.imshow(data, cmap=cmap, norm=norm, interpolation="none")
    plt.tick_params(axis="both", labelsize=TICK_SIZE)

    cbar = plt.colorbar(im)
    cbar.set_label(label, fontsize=LABEL_SIZE)
    cbar.ax.tick_params(labelsize=TICK_SIZE)


def _prefix(name: str, dir: str = FIGURES) -> str:
    return f"{dir}/{name}-" if name else ""


def _title(title: str, name: str) -> None:
    if name:
        plt.title(title + f" ({name})", fontsize=TITLE_SIZE)
    else:
        plt.title(title, fontsize=TITLE_SIZE)


def _xlabel(label: str) -> None:
    plt.xlabel(label, fontsize=LABEL_SIZE)


def _ylabel(label: str) -> None:
    plt.ylabel(label, fontsize=LABEL_SIZE)


def _save_and_close(name: str, figtype: str, dir: str = FIGURES) -> None:
    ax = plt.gca()
    _, labels = ax.get_legend_handles_labels()
    if labels:
        plt.legend(
            frameon=True,
            fancybox=False,
            shadow=True,
            fontsize=LEGEND_SIZE,
        )

    plt.savefig(_prefix(name, dir) + figtype + ".png", dpi=DPI, bbox_inches="tight")
    plt.close()


def _plot_heatmap(heatmap: ArrayType, name: str, dir: str = FIGURES) -> None:
    heatmap = np.asarray(heatmap)
    im = plt.imshow(heatmap, cmap=CMAP, norm=LogNorm(vmin=1e-4, vmax=1e0), interpolation="none")
    plt.colorbar(im, label="Occurrence Probability")

    _title("Modification Heatmap", name)
    plt.xticks(range(7), ["A", "C", "G", "U", "del", "ins", "term"], fontsize=TICK_SIZE)
    plt.yticks(range(4), ["A", "C", "G", "U"], fontsize=TICK_SIZE)

    figtype = "heatmap"
    _save_and_close(name, figtype, dir)


def _plot_read_hist(reads: ArrayType, name: str, dir: str = FIGURES) -> None:
    reads = np.asarray(reads)
    with np.errstate(divide="ignore"):
        lr = np.where(reads == 0, -1, np.log10(reads))

    hist_result = plt.hist(lr, bins=100)
    counts: np.ndarray = np.asarray(hist_result[0])
    patches = hist_result[2]
    plt.grid(axis="y", alpha=0.5)

    norm = plt.Normalize(float(counts.min()), float(counts.max()))
    for count, patch in zip(counts, patches):  # type: ignore[arg-type]
        patch.set_facecolor(plt.get_cmap("RdPu")(norm(count)))

    _title("Read distribution", name)
    _xlabel("log10 Read Depth")
    _ylabel("Count")

    figtype = "reads-hist"
    _save_and_close(name, figtype, dir)


def _plot_cumulative_reads(
    reads: ArrayType, name: str, block: int = 100, dir: str = FIGURES
) -> None:
    reads = np.asarray(reads)
    n = reads.shape[0] // block
    reads = reads[: n * block]
    reads = reads.reshape(n, block).mean(1)

    with np.errstate(divide="ignore"):
        lr = np.where(reads == 0, -1, np.log10(reads))

    _line_plot(lr)

    _title("Cumulative reads", name)
    _xlabel("Sequence index")
    _ylabel("Read Depth")

    figtype = "cumulative-reads"
    _save_and_close(name, figtype, dir)


def _plot_termination(term: ArrayType, name: str, dir: str = FIGURES) -> None:
    term = np.asarray(term)
    quot = term.sum(axis=1)[:, None]
    term = np.divide(
        term,
        quot,
        where=(quot > 0),
        out=term,
    )
    term = np.mean(term, axis=0)

    _line_plot(term)
    plt.ylim(0, 1.05)

    _title("Termination by position", name)
    _xlabel("Residue")
    _ylabel("Termination density")

    figtype = "termination"
    _save_and_close(name, figtype, dir)


def _plot_profile(
    reactivity: np.ndarray,
    error: Union[np.ndarray, None] = None,
    name: str = "",
    dir: str = FIGURES,
) -> None:
    if error is not None:
        _line_plot(reactivity, label="Reactivity")
        _line_plot(-error, plt.get_cmap("PuBu"), label="Stat. Error")
    else:
        _line_plot(reactivity)

    _title("Reactivity Profile", name)
    _xlabel("Residue")
    _ylabel("Reactivity")

    figtype = "profile"
    _save_and_close(name, figtype)


def plot_profiles(
    reactivities: list[np.ndarray],
    names: list[str],
) -> None:
    """Plot multiple reactivity profiles overlaid."""
    for ix in range(len(names)):
        _line_plot(reactivities[ix], cmap=CMAPS[ix], label=names[ix])

    name = ""

    _title("Reactivity Profile", name)
    _xlabel("Residue")
    _ylabel("Reactivity")

    figtype = "profile"
    _save_and_close(name, figtype)


def _plot_error(
    error: np.ndarray,
    name: str,
    dir: str = FIGURES,
) -> None:
    _line_plot(error)

    _title("Error Profile", name)
    _xlabel("Residue")
    _ylabel("Error")

    figtype = "error"
    _save_and_close(name, figtype, dir)


def _plot_multiple_examples(
    reactivity: np.ndarray,
    name: str,
    dir: str = FIGURES,
) -> None:
    masked = np.ma.masked_where(np.isnan(reactivity), reactivity)

    cmap = plt.get_cmap("RdPu").copy()
    cmap.set_bad(color="grey")
    im = plt.imshow(masked, cmap=cmap, vmin=0, vmax=1, interpolation="none")
    plt.colorbar(im, label="Reactivity")

    _title("Heatmap of profiles", name)
    _xlabel("Residue")
    _ylabel("Sequence Index")

    figtype = "examples"
    _save_and_close(name, figtype, dir)


def _plot_examples(
    reactivity: ArrayType,
    reads: ArrayType,
    error: ArrayType,
    name: str,
    num: int = 250,
    dir: str = FIGURES,
) -> None:
    reactivity = np.asarray(reactivity)
    error = np.asarray(error)
    if reactivity.shape[0] == 1:
        _plot_profile(reactivity[0], error[0], name, dir)
    else:
        _plot_multiple_examples(reactivity[:num], name, dir)


def _plot_coverage(coverage: ArrayType, reads: ArrayType, name: str, dir: str = FIGURES) -> None:
    coverage = np.asarray(coverage)
    reads = np.asarray(reads)
    _line_plot(coverage / reads.mean())
    plt.ylim(0, 1.05)

    _title("Coverage by position", name)
    _xlabel("Residue")
    _ylabel("Fraction of reads")

    figtype = "coverage"
    _save_and_close(name, figtype, dir)


def _plot_variance(values: np.ndarray, name: str, dir: str = FIGURES) -> None:
    _line_plot(values)

    _title("Variance", name)
    _xlabel("Residue")
    _ylabel("Variance")

    figtype = "variance"
    _save_and_close(name, figtype, dir)


def _plot_pairwise_coverage(values: np.ndarray, name: str, dir: str = FIGURES) -> None:
    vlow = max(values.min(), 1e-4)
    vhigh = 1
    norm = LogNorm(vmin=vlow, vmax=vhigh)

    _matrix_plot(values, norm, "Pairwise Correlation")

    _title("Pairwise Coverage", name)
    _xlabel("Residue")
    _ylabel("Residue")

    figtype = "pairwise-coverage"
    _save_and_close(name, figtype, dir)


def _plot_correlation(values: np.ndarray, name: str, dir: str = FIGURES) -> None:
    vhigh = 1
    norm = SymLogNorm(linthresh=1e-3 * vhigh, vmin=-vhigh, vmax=vhigh)

    _matrix_plot(values, norm, "Pairwise Correlation", cmap="PiYG")

    _title("Correlation", name)
    _xlabel("Residue")
    _ylabel("Residue")

    figtype = "correlation"
    _save_and_close(name, figtype, dir)


def _plot_mi(
    values: np.ndarray,
    name: str,
    dir: str = FIGURES,
) -> None:
    mask = ~np.isnan(values)
    vlow: float = float(np.percentile(values[mask], 95))
    vhigh: float = float(np.percentile(values[mask], 99))

    vlow = max(vlow, 1e-6)
    vhigh = max(vhigh, 1e-6)

    norm = LogNorm(vmin=vlow, vmax=vhigh)

    _matrix_plot(values, norm, "Mutual information")

    _title("Mutual Information", name)
    _xlabel("Residue")
    _ylabel("Residue")

    figtype = "mutual-information"
    _save_and_close(name, figtype, dir)


def plot_snr_scaling(
    mod: ProbingData,
    nomod: "Union[ProbingData, None]",
    combined: ProbingData,
    name: str,
    dir: str = FIGURES,
) -> None:
    """Plot expected mean SNR as a function of relative total read depth.

    When nomod is present, three curves show how extra reads could be allocated:
    - Modified: all extra reads go to the modified condition
    - Unmodified: all extra reads go to the unmodified condition
    - Both: extra reads split proportionally between conditions
    """
    os.makedirs(dir, exist_ok=True)

    x = np.geomspace(1e-4, 10, 200)
    reactivity = np.asarray(combined.reactivity)
    mod_err = np.asarray(mod.error)
    mod_reads = float(np.asarray(mod.reads).max())

    if nomod is not None:
        nomod_err = np.asarray(nomod.error)
        nomod_reads = float(np.asarray(nomod.reads).max())
        total_reads = mod_reads + nomod_reads
        f_mod = mod_reads / total_reads
        f_nomod = nomod_reads / total_reads

    def _mean_snr(mod_scale: float, nomod_scale: float) -> float:
        se2 = mod_err**2 / mod_scale + (nomod_err**2 / nomod_scale if nomod is not None else 0)
        se = np.sqrt(se2)
        snr = np.divide(reactivity, se, where=(se > 0), out=np.zeros_like(reactivity))
        return float(np.nanmean(snr, axis=-1).mean(axis=0))

    k = x
    if nomod is not None:
        # Both: scale both proportionally
        snr_both = np.array([_mean_snr(ki, ki) for ki in k])
        plt.plot(k, snr_both, color="grey", linewidth=2, linestyle="--", label="Both")

        # Optimal: maximize mean SNR over allocation fraction
        from scipy.optimize import minimize_scalar

        def _optimal_snr(ki: float) -> float:
            def neg_snr(f: float) -> float:
                n_m = ki * total_reads * f
                n_n = ki * total_reads * (1 - f)
                return -_mean_snr(n_m / mod_reads, n_n / nomod_reads)
            res = minimize_scalar(neg_snr, bounds=(1e-3, 1 - 1e-3), method="bounded")
            return -res.fun

        snr_opt = np.array([_optimal_snr(ki) for ki in k])
        plt.plot(k, snr_opt, color="red", linewidth=2, label="Optimal")

        # Modified: scale mod only, remap x to relative total
        snr_mod = np.array([_mean_snr(ki, 1.0) for ki in k])
        x_mod = k * f_mod + f_nomod
        cmap_mod = plt.get_cmap("RdPu")
        plt.plot(x_mod, snr_mod, color=cmap_mod(0.7), linewidth=2, label="Modified")

        # Unmodified: scale nomod only, remap x to relative total
        snr_nomod = np.array([_mean_snr(1.0, ki) for ki in k])
        x_nomod = f_mod + k * f_nomod
        cmap_nomod = plt.get_cmap("PuBu")
        plt.plot(x_nomod, snr_nomod, color=cmap_nomod(0.7), linewidth=2, label="Unmodified")
    else:
        snr_mod = np.array([_mean_snr(ki, 1.0) for ki in k])
        cmap_mod = plt.get_cmap("RdPu")
        plt.plot(k, snr_mod, color=cmap_mod(0.7), linewidth=2, label="Modified")

    plt.axvline(1.0, color="grey", linestyle="--", linewidth=1, alpha=0.7)
    plt.xscale("log")
    plt.xlim(0.1, 10)
    plt.grid(axis="y", alpha=0.5)
    plt.tick_params(axis="both", labelsize=TICK_SIZE)

    _title("SNR vs Read Depth", name)
    _xlabel("Relative Total Read Depth")
    _ylabel("Mean SNR")

    _save_and_close(name, "snr-scaling", dir)


def plot_all(
    data: ProbingData,
    name: str,
    dir: str = FIGURES,
) -> None:
    """Generate all standard plots for probing data."""
    if not isinstance(data, ProbingData):
        raise TypeError(f"Expected ProbingData, got {type(data).__name__}")

    os.makedirs(dir, exist_ok=True)

    _plot_heatmap(data.heatmap, name, dir)
    _plot_examples(data.reactivity, data.reads, data.error, name, dir=dir)
    _plot_termination(data.terminations, name, dir)
    _plot_coverage(data.coverage, data.reads, name, dir)

    if not data.single():
        _plot_read_hist(data.reads, name, dir)
        _plot_cumulative_reads(data.reads, name, dir=dir)

    if data.mi is not None:
        for ix in range(data.size()):
            _plot_mi(data.mi[ix], name, dir)

    if data.covariance is not None:
        for ix in range(data.size()):
            _plot_correlation(data.covariance[ix], name, dir)


def main():
    """CLI entry point for cmuts-plot."""
    import argparse

    import h5py

    parser = argparse.ArgumentParser(
        prog="cmuts-plot", description="Generate plots from reactivity data"
    )
    parser.add_argument("file", help="HDF5 file with reactivity data")
    parser.add_argument("--group", default="", help="Group name in HDF5 file")
    parser.add_argument("-o", "--out", default=FIGURES, help="Output directory")
    args = parser.parse_args()

    if not os.path.exists(args.file):
        raise FileNotFoundError(f"File not found: {args.file}")

    with h5py.File(args.file, "r") as f:
        data = ProbingData.load(args.group, f)

    plot_all(data, args.group, args.out)
    print(f"Plots saved to {args.out}/")


if __name__ == "__main__":
    main()
