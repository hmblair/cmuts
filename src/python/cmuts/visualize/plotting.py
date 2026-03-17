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

    When nomod is present, three curves are shown:
    - Modified: all extra reads go to the modified condition
    - Unmodified: all extra reads go to the unmodified condition
    - Pareto: optimal allocation of reads between conditions

    When only a modified condition is present, a single curve is shown.
    """
    from matplotlib.transforms import ScaledTranslation

    os.makedirs(dir, exist_ok=True)

    reactivity = np.asarray(combined.reactivity)
    mod_err = np.asarray(mod.error)
    mod_reads = float(np.asarray(mod.reads).max())

    prior = 0.001
    mod_err2 = np.maximum(mod_err**2, prior * (1 - prior) / mod_reads)

    if nomod is not None:
        nomod_err = np.asarray(nomod.error)
        nomod_reads = float(np.asarray(nomod.reads).max())
        total_reads = mod_reads + nomod_reads
        nomod_err2 = np.maximum(nomod_err**2, prior * (1 - prior) / nomod_reads)
    else:
        nomod_err2 = None
        total_reads = mod_reads

    # SEM of mean SNR at current depth. Error bands at projected depths
    # scale this proportionally by curve(k)/curve(1).
    current_se2 = mod_err2.copy()
    if nomod_err2 is not None:
        current_se2 = current_se2 + nomod_err2
    current_se = np.sqrt(current_se2)
    current_snr = np.where(current_se > 0, reactivity / current_se, 0.0)
    _flat = current_snr.ravel()
    _valid = _flat[~np.isnan(_flat)]
    snr_at_1 = _valid.mean() if len(_valid) > 0 else 1.0
    snr_sem_at_1 = _valid.std() / np.sqrt(len(_valid)) if len(_valid) > 1 else 0.0

    def _mean_snr_vec(mod_scales: np.ndarray, nomod_scales: np.ndarray) -> np.ndarray:
        """Compute mean SNR for an array of scale factors. Returns shape (N,)."""
        se2 = mod_err2[None] / mod_scales[:, None, None]
        if nomod_err2 is not None:
            se2 = se2 + nomod_err2[None] / nomod_scales[:, None, None]
        se = np.sqrt(se2)
        snr = np.where(se > 0, reactivity[None] / se, 0.0)
        return np.nanmean(snr, axis=-1).mean(axis=-1)

    def _snr_band(snr_curve: np.ndarray) -> np.ndarray:
        """Compute ±1 SEM band by scaling the current SEM proportionally."""
        ratio = np.where(snr_at_1 > 0, snr_curve / snr_at_1, 0.0)
        return np.abs(ratio) * snr_sem_at_1

    def _trim_leading_zeros(snr: np.ndarray) -> slice:
        """Trim leading near-zero region, keeping one point before the visible rise."""
        nz = np.nonzero(snr > 0.01)[0]
        start = max(nz[0] - 1, 0) if len(nz) > 0 else 0
        return slice(start, None)

    xi = np.geomspace(0.1, 10, 10000)

    if nomod is not None:
        # Modified: extra reads go to mod, nomod stays at current depth
        mod_scales = np.maximum((xi * total_reads - nomod_reads) / mod_reads, 1e-10)
        snr_mod = _mean_snr_vec(mod_scales, np.ones_like(mod_scales))
        sem_mod = _snr_band(snr_mod)
        s = _trim_leading_zeros(snr_mod)
        mod_color = plt.get_cmap("RdPu")(0.7)
        plt.fill_between(
            xi[s],
            (snr_mod - sem_mod)[s],
            (snr_mod + sem_mod)[s],
            color=mod_color,
            alpha=0.2,
        )
        plt.plot(xi[s], snr_mod[s], color=mod_color, linewidth=2, label="Modified")

        # Unmodified: extra reads go to nomod, mod stays at current depth
        nomod_scales = np.maximum((xi * total_reads - mod_reads) / nomod_reads, 1e-10)
        snr_nomod = _mean_snr_vec(np.ones_like(nomod_scales), nomod_scales)
        sem_nomod = _snr_band(snr_nomod)
        s = _trim_leading_zeros(snr_nomod)
        nomod_color = plt.get_cmap("PuBu")(0.7)
        plt.fill_between(
            xi[s],
            (snr_nomod - sem_nomod)[s],
            (snr_nomod + sem_nomod)[s],
            color=nomod_color,
            alpha=0.2,
        )
        plt.plot(xi[s], snr_nomod[s], color=nomod_color, linewidth=2, label="Unmodified")

        # Set ylim based on mod/nomod curves before Pareto expands it
        curve_max: float = max(np.nanmax(snr_mod + sem_mod), np.nanmax(snr_nomod + sem_nomod))
        plt.ylim(0, curve_max * 1.1)

        # Pareto: best allocation at each relative total depth
        fracs = np.linspace(0.01, 0.99, 200)
        snr_pareto = np.empty(len(xi))
        chunk = 500
        for i in range(0, len(xi), chunk):
            xi_c = xi[i : i + chunk]
            ms = xi_c[:, None] * total_reads * fracs[None, :] / mod_reads
            ns = xi_c[:, None] * total_reads * (1 - fracs[None, :]) / nomod_reads
            snr_grid = np.column_stack(
                [_mean_snr_vec(ms[:, j], ns[:, j]) for j in range(len(fracs))]
            )
            snr_pareto[i : i + chunk] = snr_grid.max(axis=1)

        sem_pareto = _snr_band(snr_pareto)
        pareto_upper = snr_pareto + sem_pareto

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
        snr_mod = _mean_snr_vec(xi, np.ones_like(xi))
        sem_mod = _snr_band(snr_mod)
        mod_color = plt.get_cmap("RdPu")(0.7)
        plt.fill_between(
            xi,
            snr_mod - sem_mod,
            snr_mod + sem_mod,
            color=mod_color,
            alpha=0.2,
        )
        plt.plot(xi, snr_mod, color=mod_color, linewidth=2, label="Modified")

    plt.axvline(1.0, color="grey", linestyle="--", linewidth=1, alpha=0.7)
    plt.xscale("log")
    plt.xlim(0.1, 10)
    plt.grid(axis="y", alpha=0.5)
    plt.tick_params(axis="both", labelsize=TICK_SIZE)

    # Hatch infeasible region above Pareto curve
    if nomod is not None:
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
