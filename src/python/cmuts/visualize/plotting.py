import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, SymLogNorm
from cmuts.internal import ProbingData
from typing import Union


# Figures directory


FIGURES = "figures"


# Font and color for plotting


plt.rcParams['font.family'] = 'Helvetica'
LABEL_SIZE = 14
TITLE_SIZE = 15
TICK_SIZE = 13
LEGEND_SIZE = 12
DPI = 400

CMAP = "RdPu"
CMAPS = [
    plt.cm.RdPu,
    plt.cm.PiYG,
]


def _line_plot(data: np.ndarray, cmap = plt.cm.RdPu, label: str = "") -> None:

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
    plt.tick_params(axis='both', labelsize=TICK_SIZE)
    plt.xlim(0, n)


def _matrix_plot(data: np.ndarray, norm, label: str, cmap: str = CMAP) -> None:

    im = plt.imshow(data, cmap=cmap, norm=norm, interpolation='none')
    plt.tick_params(axis='both', labelsize=TICK_SIZE)

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


def _plot_heatmap(heatmap: np.ndarray, name: str, dir: str = FIGURES) -> None:

    im = plt.imshow(heatmap, cmap=CMAP, norm=LogNorm(vmin=1E-4, vmax=1E0), interpolation='none')
    plt.colorbar(im, label="Occurrence Probability")

    _title("Modification Heatmap", name)
    plt.xticks(range(7), ["A", "C", "G", "U", "del", "ins", "term"], fontsize=TICK_SIZE)
    plt.yticks(range(4), ["A", "C", "G", "U"], fontsize=TICK_SIZE)

    figtype = "heatmap"
    _save_and_close(name, figtype, dir)


def _plot_read_hist(reads: np.ndarray, name: str, dir: str = FIGURES) -> None:

    with np.errstate(divide='ignore'):
        lr = np.where(reads == 0, -1, np.log10(reads))

    counts, bins, patches = plt.hist(lr, bins=100)
    plt.grid(axis="y", alpha=0.5)

    norm = plt.Normalize(counts.min(), counts.max())
    for count, patch in zip(counts, patches):
        patch.set_facecolor(plt.cm.RdPu(norm(count)))

    _title("Read distribution", name)
    _xlabel("log10 Read Depth")
    _ylabel("Count")

    figtype = "reads-hist"
    _save_and_close(name, figtype, dir)


def _plot_cumulative_reads(reads: np.ndarray, name: str, block: int = 100, dir: str = FIGURES) -> None:

    n = reads.shape[0] // block
    reads = reads[:n * block]
    reads = reads.reshape(n, block).mean(1)

    with np.errstate(divide='ignore'):
        lr = np.where(reads == 0, -1, np.log10(reads))

    _line_plot(lr)

    _title("Cumulative reads", name)
    _xlabel("Sequence index")
    _ylabel("Read Depth")

    figtype = "cumulative-reads"
    _save_and_close(name, figtype, dir)


def _plot_termination(term: np.ndarray, name: str, dir: str = FIGURES) -> None:

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
        _line_plot(-error, plt.cm.PuBu, label="Stat. Error")
    else:
        _line_plot(reactivity)

    _title("Reactivity Profile", name)
    _xlabel("Residue")
    _ylabel("Reactivity")

    figtype = "profile"
    _save_and_close(name, figtype)


def _plot_profiles(
    reactivities: list[np.ndarray],
    names: list[str],
) -> None:

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

    cmap = plt.cm.RdPu.copy()
    cmap.set_bad(color='grey')
    im = plt.imshow(masked, cmap=cmap, vmin=0, vmax=1, interpolation='none')
    plt.colorbar(im, label="Reactivity")

    _title("Heatmap of profiles", name)
    _xlabel("Residue")
    _ylabel("Sequence Index")

    figtype = "examples"
    _save_and_close(name, figtype, dir)


def _plot_examples(
    reactivity: np.ndarray,
    reads: np.ndarray,
    error: np.ndarray,
    name: str,
    num: int = 250,
    dir: str = FIGURES,
) -> None:

    if reactivity.shape[0] == 1:
        _plot_profile(reactivity[0], error[0], name, dir)
    else:
        _plot_multiple_examples(reactivity[:num], name, dir)


def _plot_coverage(coverage: np.ndarray, reads: np.ndarray, name: str, dir: str = FIGURES) -> None:

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

    vlow = values.min()
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
    norm = SymLogNorm(linthresh=1E-3*vhigh, vmin=-vhigh, vmax=vhigh)

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
    vlow = np.percentile(values[mask], 95).astype(float)
    vhigh = np.percentile(values[mask], 99).astype(float)
    vlow = max(vlow, 1E-6)
    norm = LogNorm(vmin=vlow, vmax=vhigh)

    _matrix_plot(values, norm, "Mutual information")

    _title("Mutual Information", name)
    _xlabel("Residue")
    _ylabel("Residue")

    figtype = "mutual-information"
    _save_and_close(name, figtype, dir)


def all(
    data: ProbingData,
    name: str,
    dir: str = FIGURES,
) -> None:

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
