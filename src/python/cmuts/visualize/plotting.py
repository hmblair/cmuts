import os
from typing import Any, Union

import dask.array as da
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm, SymLogNorm

from cmuts.internal import ProbingData, SNRCurves

from . import _transforms

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
    if name:
        return os.path.join(dir, f"{name}-")
    return os.path.join(dir, "")


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
    os.makedirs(dir, exist_ok=True)
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
    lr = _transforms.read_log_depth(reads)

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


def _plot_reads_per_block(
    reads: ArrayType, name: str, nblocks: int = 100, dir: str = FIGURES
) -> None:
    """Plot mean reads per reference, references split into a fixed
    number of equal-sized blocks in their FASTA order.
    """
    means = _transforms.reads_per_block(reads, nblocks)

    _line_plot(means)
    _title("Mean reads per reference bin", name)
    _xlabel("Reference bin")
    _ylabel("Mean reads")

    figtype = "reads-per-block"
    _save_and_close(name, figtype, dir)


def _plot_termination(term: ArrayType, name: str, dir: str = FIGURES) -> None:
    term = _transforms.termination_density(term)

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
    _save_and_close(name, figtype, dir)


def plot_profiles(
    reactivities: list[np.ndarray],
    names: list[str],
    dir: str = FIGURES,
) -> None:
    """Plot multiple reactivity profiles overlaid."""
    for ix in range(len(names)):
        _line_plot(reactivities[ix], cmap=CMAPS[ix], label=names[ix])

    name = ""

    _title("Reactivity Profile", name)
    _xlabel("Residue")
    _ylabel("Reactivity")

    figtype = "profile"
    _save_and_close(name, figtype, dir)


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
    _line_plot(_transforms.coverage_fraction(coverage, reads))
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
    vlow, vhigh = _transforms.pairwise_coverage_bounds(values)
    norm = LogNorm(vmin=vlow, vmax=vhigh)

    _matrix_plot(values, norm, "Pairwise Correlation")

    _title("Pairwise Coverage", name)
    _xlabel("Residue")
    _ylabel("Residue")

    figtype = "pairwise-coverage"
    _save_and_close(name, figtype, dir)


def _plot_correlation(values: np.ndarray, name: str, dir: str = FIGURES) -> None:
    norm = SymLogNorm(linthresh=_transforms.CORRELATION_LINTHRESH, vmin=-1, vmax=1)

    _matrix_plot(values, norm, "Pairwise Correlation", cmap="PiYG")

    _title("Correlation", name)
    _xlabel("Residue")
    _ylabel("Residue")

    figtype = "correlation"
    _save_and_close(name, figtype, dir)


def _indexed_name(name: str, ix: int, total: int) -> str:
    if total <= 1:
        return name
    suffix = f"{ix + 1}"
    return f"{name}-{suffix}" if name else suffix


def _plot_mi(
    values: np.ndarray,
    name: str,
    dir: str = FIGURES,
) -> None:
    vlow, vhigh = _transforms.mi_bounds(values)

    norm = LogNorm(vmin=vlow, vmax=vhigh)

    _matrix_plot(values, norm, "Mutual information")

    _title("Mutual Information", name)
    _xlabel("Residue")
    _ylabel("Residue")

    figtype = "mutual-information"
    _save_and_close(name, figtype, dir)


def plot_snr_scaling(curves: SNRCurves, name: str, dir: str = FIGURES) -> None:
    """Render precomputed SNR-vs-read-depth curves (see cmuts.compute_snr_curves)."""
    from matplotlib.transforms import ScaledTranslation

    os.makedirs(dir, exist_ok=True)

    def _trim(snr: np.ndarray) -> slice:
        """Drop the leading near-zero region, keeping one point before the rise."""
        nz = np.nonzero(snr > 0.01)[0]
        start = max(nz[0] - 1, 0) if len(nz) > 0 else 0
        return slice(start, None)

    xi = curves.xi
    mod_color = plt.get_cmap("RdPu")(0.7)

    if curves.nomod is not None:
        # nomod_sem, pareto and pareto_sem are populated together with nomod.
        assert curves.nomod_sem is not None
        assert curves.pareto is not None and curves.pareto_sem is not None
        s = _trim(curves.mod)
        plt.fill_between(
            xi[s],
            (curves.mod - curves.mod_sem)[s],
            (curves.mod + curves.mod_sem)[s],
            color=mod_color,
            alpha=0.2,
        )
        plt.plot(xi[s], curves.mod[s], color=mod_color, linewidth=2, label="Modified")

        nomod_color = plt.get_cmap("PuBu")(0.7)
        s = _trim(curves.nomod)
        plt.fill_between(
            xi[s],
            (curves.nomod - curves.nomod_sem)[s],
            (curves.nomod + curves.nomod_sem)[s],
            color=nomod_color,
            alpha=0.2,
        )
        plt.plot(xi[s], curves.nomod[s], color=nomod_color, linewidth=2, label="Unmodified")

        curve_max: float = max(
            np.nanmax(curves.mod + curves.mod_sem), np.nanmax(curves.nomod + curves.nomod_sem)
        )
        plt.ylim(0, curve_max * 1.1)

        pareto_upper = curves.pareto + curves.pareto_sem
        dy = ScaledTranslation(0, 2.0 / 72, plt.gcf().dpi_scale_trans)
        pareto_transform = plt.gca().transData + dy
        plt.plot(
            xi,
            pareto_upper,
            color="black",
            linewidth=1,
            linestyle=(0, (3, 2)),
            label="Pareto",
            zorder=6,
            transform=pareto_transform,
        )
    else:
        plt.fill_between(
            xi,
            curves.mod - curves.mod_sem,
            curves.mod + curves.mod_sem,
            color=mod_color,
            alpha=0.2,
        )
        plt.plot(xi, curves.mod, color=mod_color, linewidth=2, label="Modified")

    plt.axvline(1.0, color="grey", linestyle="--", linewidth=1, alpha=0.7)
    plt.xscale("log")
    plt.xlim(xi[0], xi[-1])
    plt.grid(axis="y", alpha=0.5)
    plt.tick_params(axis="both", labelsize=TICK_SIZE)

    if curves.nomod is not None:
        ax = plt.gca()
        ylim = ax.get_ylim()
        ax.fill_between(
            xi,
            pareto_upper,
            ylim[1] * 2,
            color="none",
            edgecolor="silver",
            linewidth=0,
            hatch="/////",
            alpha=1.0,
            zorder=5,
            transform=pareto_transform,
        )
        for spine in ax.spines.values():
            spine.set_zorder(10)
        ax.tick_params(zorder=10)
        ax.set_ylim(ylim)

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
        _plot_reads_per_block(data.reads, name, dir=dir)

    if data.mi is not None:
        for ix in range(data.size()):
            _plot_mi(data.mi[ix], _indexed_name(name, ix, data.size()), dir)

    if data.covariance is not None:
        for ix in range(data.size()):
            _plot_correlation(data.covariance[ix], _indexed_name(name, ix, data.size()), dir)


def main():
    """CLI entry point for cmuts-plot: render every figure from a reactivity h5."""
    import argparse

    import h5py

    parser = argparse.ArgumentParser(
        prog="cmuts-plot", description="Generate plots from a reactivity HDF5 file"
    )
    parser.add_argument("file", help="HDF5 file written by cmuts normalize")
    parser.add_argument("--group", default=None, help="Plot only this group (default: all)")
    parser.add_argument("-o", "--out", default=FIGURES, help="Output directory")
    args = parser.parse_args()

    if not os.path.exists(args.file):
        raise FileNotFoundError(f"File not found: {args.file}")

    with h5py.File(args.file, "r") as f:
        if args.group is not None:
            groups = [args.group]
        else:
            groups = [k for k in f if isinstance(f[k], h5py.Group)] or [""]
        for group in groups:
            plot_all(ProbingData.load(group, f), group, args.out)
            curves = SNRCurves.load(group, f)
            if curves is not None:
                plot_snr_scaling(curves, group, args.out)
    print(f"Plots saved to {args.out}/")


if __name__ == "__main__":
    main()
