"""Plotly implementations of cmuts plots.

Each function returns a ``go.Figure`` for interactive display (e.g. in Gradio).
The matplotlib module (``plotting.py``) remains the CLI default and saves PNGs.
"""

from typing import Any, Optional

import numpy as np
import plotly.colors as pc
import plotly.graph_objects as go

from cmuts.internal import ProbingData

# Styling constants matching the matplotlib versions
FONT_FAMILY = "Helvetica"
LABEL_SIZE = 14
TITLE_SIZE = 15
TICK_SIZE = 13
LEGEND_SIZE = 12

# Pre-sampled colors from colorscales
_RDPU_FILL = pc.sample_colorscale("RdPu", [0.3])[0]
_RDPU_LINE = pc.sample_colorscale("RdPu", [0.8])[0]
_RDPU_070 = pc.sample_colorscale("RdPu", [0.7])[0]
_PUBU_FILL = pc.sample_colorscale("PuBu", [0.3])[0]
_PUBU_LINE = pc.sample_colorscale("PuBu", [0.8])[0]
_PUBU_070 = pc.sample_colorscale("PuBu", [0.7])[0]


def _rgba(rgb_str: str, alpha: float) -> str:
    return rgb_str.replace("rgb(", "rgba(").replace(")", f", {alpha})")


def _base_layout(title: str, xlabel: str = "", ylabel: str = "") -> dict[str, Any]:
    return {
        "font": {"family": FONT_FAMILY},
        "title": {"text": title, "font": {"size": TITLE_SIZE}},
        "xaxis": {
            "title": {"text": xlabel, "font": {"size": LABEL_SIZE}},
            "tickfont": {"size": TICK_SIZE},
        },
        "yaxis": {
            "title": {"text": ylabel, "font": {"size": LABEL_SIZE}},
            "tickfont": {"size": TICK_SIZE},
            "gridcolor": "rgba(0,0,0,0.15)",
        },
        "plot_bgcolor": "white",
        "legend": {"font": {"size": LEGEND_SIZE}},
    }


def _line_trace(
    fig: go.Figure,
    data: np.ndarray,
    fill_color: str = _RDPU_FILL,
    line_color: str = _RDPU_LINE,
    label: str = "",
) -> None:
    data = np.asarray(data)
    x = np.arange(data.shape[0])
    fig.add_trace(
        go.Scatter(
            x=x,
            y=data,
            fill="tozeroy",
            fillcolor=_rgba(fill_color, 0.5),
            line={"color": line_color, "width": 1},
            name=label,
            showlegend=bool(label),
        )
    )


def _symlog(x: np.ndarray, linthresh: float) -> np.ndarray:
    return np.sign(x) * np.log10(1.0 + np.abs(x) / linthresh)


# ---------------------------------------------------------------------------
# Public plot functions — each returns a go.Figure
# ---------------------------------------------------------------------------


_HEATMAP_NTS = ["A", "C", "G", "U"]
_HEATMAP_MODS = ["A", "C", "G", "U", "del", "ins", "term"]


def plot_heatmap(heatmap: np.ndarray, name: str = "") -> go.Figure:
    heatmap = np.asarray(heatmap)
    with np.errstate(divide="ignore"):
        log_heatmap = np.where(heatmap > 0, np.log10(heatmap), np.nan)

    # Build descriptive hover text per cell
    hover_text = []
    for i, nt in enumerate(_HEATMAP_NTS):
        row = []
        for j, mod in enumerate(_HEATMAP_MODS):
            val = heatmap[i, j]
            prob = f"{val:.4e}" if val > 0 else "0"
            if mod in _HEATMAP_NTS and mod == nt:
                row.append(f"Match ({nt})<br>Probability: {prob}")
            elif mod in _HEATMAP_NTS:
                row.append(f"Mismatch {nt} \u2192 {mod}<br>Probability: {prob}")
            elif mod == "del":
                row.append(f"Deletion of {nt}<br>Probability: {prob}")
            elif mod == "ins":
                row.append(f"Insertion at {nt}<br>Probability: {prob}")
            else:
                row.append(f"Termination at {nt}<br>Probability: {prob}")
        hover_text.append(row)

    title = f"Modification Heatmap ({name})" if name else "Modification Heatmap"
    fig = go.Figure(
        data=go.Heatmap(
            z=log_heatmap,
            x=_HEATMAP_MODS,
            y=_HEATMAP_NTS,
            colorscale="RdPu",
            zmin=-4,
            zmax=0,
            text=hover_text,
            hoverinfo="text",
            colorbar={
                "title": "Probability",
                "tickvals": [-4, -3, -2, -1, 0],
                "ticktext": [
                    "10\u207b\u2074",
                    "10\u207b\u00b3",
                    "10\u207b\u00b2",
                    "10\u207b\u00b9",
                    "10\u2070",
                ],
            },
        )
    )

    # Cell outlines
    for i in range(len(_HEATMAP_NTS)):
        for j in range(len(_HEATMAP_MODS)):
            fig.add_shape(
                type="rect",
                x0=j - 0.5,
                x1=j + 0.5,
                y0=i - 0.5,
                y1=i + 0.5,
                line={"color": "black", "width": 1},
                layer="above",
            )

    fig.update_layout(
        font={"family": FONT_FAMILY},
        title={"text": title, "font": {"size": TITLE_SIZE}},
        xaxis_title="Modification Type",
        xaxis={"showgrid": False, "zeroline": False, "constrain": "domain"},
        yaxis_title="Reference Nucleotide",
        yaxis={
            "autorange": "reversed",
            "showgrid": False,
            "zeroline": False,
            "ticklabelstandoff": 10,
            "scaleanchor": "x",
            "constrain": "domain",
        },
        template="plotly_white",
        height=300,
        width=550,
        margin={"l": 50, "r": 20, "t": 40, "b": 40},
    )
    return fig


def plot_read_hist(reads: np.ndarray, name: str = "") -> go.Figure:
    reads = np.asarray(reads)
    with np.errstate(divide="ignore"):
        lr = np.where(reads == 0, -1, np.log10(reads))

    counts, bin_edges = np.histogram(lr, bins=100)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    cmin, cmax = float(counts.min()), float(counts.max())
    if cmax == cmin:
        normed = np.zeros_like(counts, dtype=float)
    else:
        normed = (counts - cmin) / (cmax - cmin)
    colors = [pc.sample_colorscale("RdPu", [float(v)])[0] for v in normed]

    title = f"Read distribution ({name})" if name else "Read distribution"
    fig = go.Figure(
        data=go.Bar(
            x=bin_centers,
            y=counts,
            marker_color=colors,
            showlegend=False,
        )
    )
    fig.update_layout(**_base_layout(title, "log10 Read Depth", "Count"))
    fig.update_layout(bargap=0)
    return fig


def plot_cumulative_reads(reads: np.ndarray, name: str = "", block: int = 100) -> go.Figure:
    reads = np.asarray(reads)
    block = min(block, reads.shape[0])
    n = max(reads.shape[0] // block, 1)
    reads = reads[: n * block]
    reads = reads.reshape(n, block).mean(1)

    with np.errstate(divide="ignore"):
        lr = np.where(reads == 0, -1, np.log10(reads))

    title = f"Cumulative reads ({name})" if name else "Cumulative reads"
    fig = go.Figure()
    _line_trace(fig, lr)
    fig.update_layout(**_base_layout(title, "Sequence index", "Read Depth"))
    fig.update_xaxes(range=[0, len(lr)])
    return fig


def plot_termination(term: np.ndarray, name: str = "") -> go.Figure:
    term = np.asarray(term).copy()
    quot = term.sum(axis=1)[:, None]
    term = np.divide(term, quot, where=(quot > 0), out=term)
    term = np.mean(term, axis=0)

    title = f"Termination by position ({name})" if name else "Termination by position"
    fig = go.Figure()
    _line_trace(fig, term)
    fig.update_layout(**_base_layout(title, "Residue", "Termination density"))
    fig.update_xaxes(range=[0, len(term)])
    fig.update_yaxes(range=[0, 1.05])
    return fig


def plot_profile(
    reactivity: np.ndarray,
    error: Optional[np.ndarray] = None,
    name: str = "",
) -> go.Figure:
    title = f"Reactivity Profile ({name})" if name else "Reactivity Profile"
    fig = go.Figure()

    if error is not None:
        _line_trace(fig, reactivity, label="Reactivity")
        _line_trace(fig, -error, _PUBU_FILL, _PUBU_LINE, label="Stat. Error")
    else:
        _line_trace(fig, reactivity)

    fig.update_layout(**_base_layout(title, "Residue", "Reactivity"))
    fig.update_xaxes(range=[0, len(reactivity)])
    return fig


def plot_profiles(
    reactivities: list[np.ndarray],
    names: list[str],
) -> go.Figure:
    color_sets = [
        (_RDPU_FILL, _RDPU_LINE),
        (pc.sample_colorscale("PiYG", [0.3])[0], pc.sample_colorscale("PiYG", [0.8])[0]),
    ]
    fig = go.Figure()
    for ix, (reac, label) in enumerate(zip(reactivities, names)):
        fill_c, line_c = color_sets[ix % len(color_sets)]
        _line_trace(fig, reac, fill_c, line_c, label=label)

    fig.update_layout(**_base_layout("Reactivity Profile", "Residue", "Reactivity"))
    max_len = max(len(r) for r in reactivities)
    fig.update_xaxes(range=[0, max_len])
    return fig


def plot_error(error: np.ndarray, name: str = "") -> go.Figure:
    title = f"Error Profile ({name})" if name else "Error Profile"
    fig = go.Figure()
    _line_trace(fig, error)
    fig.update_layout(**_base_layout(title, "Residue", "Error"))
    fig.update_xaxes(range=[0, len(error)])
    return fig


def plot_multiple_examples(reactivity: np.ndarray, name: str = "") -> go.Figure:
    title = f"Heatmap of profiles ({name})" if name else "Heatmap of profiles"

    sentinel = -0.01
    data = reactivity.copy()
    data[np.isnan(data)] = sentinel

    eps = 0.001
    colorscale = [
        [0.0, "grey"],
        [eps, "grey"],
        [eps + 0.001, pc.sample_colorscale("RdPu", [0.0])[0]],
        [1.0, pc.sample_colorscale("RdPu", [1.0])[0]],
    ]

    fig = go.Figure(
        data=go.Heatmap(
            z=data,
            colorscale=colorscale,
            zmin=sentinel,
            zmax=1,
            colorbar={"title": "Reactivity"},
        )
    )
    fig.update_layout(**_base_layout(title, "Residue", "Sequence Index"))
    fig.update_layout(yaxis_autorange="reversed")
    return fig


def plot_examples(
    reactivity: np.ndarray,
    error: np.ndarray,
    name: str = "",
    num: int = 250,
) -> go.Figure:
    reactivity = np.asarray(reactivity)
    error = np.asarray(error)
    if reactivity.shape[0] == 1:
        return plot_profile(reactivity[0], error[0], name)
    else:
        return plot_multiple_examples(reactivity[:num], name)


def plot_coverage(coverage: np.ndarray, reads: np.ndarray, name: str = "") -> go.Figure:
    coverage = np.asarray(coverage)
    reads = np.asarray(reads)
    data = coverage / reads.mean()

    title = f"Coverage by position ({name})" if name else "Coverage by position"
    fig = go.Figure()
    _line_trace(fig, data)
    fig.update_layout(**_base_layout(title, "Residue", "Fraction of reads"))
    fig.update_xaxes(range=[0, len(data)])
    fig.update_yaxes(range=[0, 1.05])
    return fig


def plot_variance(values: np.ndarray, name: str = "") -> go.Figure:
    title = f"Variance ({name})" if name else "Variance"
    fig = go.Figure()
    _line_trace(fig, values)
    fig.update_layout(**_base_layout(title, "Residue", "Variance"))
    fig.update_xaxes(range=[0, len(values)])
    return fig


def plot_pairwise_coverage(values: np.ndarray, name: str = "") -> go.Figure:
    vlow = max(float(values.min()), 1e-4)
    with np.errstate(divide="ignore"):
        log_values = np.log10(np.clip(values, vlow, 1.0))

    title = f"Pairwise Coverage ({name})" if name else "Pairwise Coverage"
    fig = go.Figure(
        data=go.Heatmap(
            z=log_values,
            colorscale="RdPu",
            zmin=np.log10(vlow),
            zmax=0,
            colorbar={"title": "Pairwise Correlation"},
        )
    )
    fig.update_layout(**_base_layout(title, "Residue", "Residue"))
    fig.update_layout(yaxis_autorange="reversed")
    return fig


def plot_correlation(values: np.ndarray, name: str = "") -> go.Figure:
    vhigh = 1.0
    linthresh = 1e-3 * vhigh
    transformed = _symlog(values, linthresh)
    t_max = _symlog(np.array([vhigh]), linthresh)[0]

    title = f"Correlation ({name})" if name else "Correlation"
    fig = go.Figure(
        data=go.Heatmap(
            z=transformed,
            colorscale="PiYG",
            zmin=-t_max,
            zmax=t_max,
            colorbar={
                "title": "Pairwise Correlation",
                "tickvals": [
                    _symlog(np.array([v]), linthresh)[0] for v in [-1, -0.1, -0.01, 0, 0.01, 0.1, 1]
                ],
                "ticktext": ["-1", "-0.1", "-0.01", "0", "0.01", "0.1", "1"],
            },
        )
    )
    fig.update_layout(**_base_layout(title, "Residue", "Residue"))
    fig.update_layout(yaxis_autorange="reversed")
    return fig


def plot_mi(values: np.ndarray, name: str = "") -> go.Figure:
    mask = ~np.isnan(values)
    vlow: float = max(float(np.percentile(values[mask], 95)), 1e-6)
    vhigh: float = max(float(np.percentile(values[mask], 99)), 1e-6)

    with np.errstate(divide="ignore"):
        log_values = np.log10(np.clip(values, vlow, vhigh))

    title = f"Mutual Information ({name})" if name else "Mutual Information"
    fig = go.Figure(
        data=go.Heatmap(
            z=log_values,
            colorscale="RdPu",
            zmin=np.log10(vlow),
            zmax=np.log10(vhigh),
            colorbar={"title": "Mutual information"},
        )
    )
    fig.update_layout(**_base_layout(title, "Residue", "Residue"))
    fig.update_layout(yaxis_autorange="reversed")
    return fig


def plot_snr_scaling(
    mod: ProbingData,
    nomod: Optional[ProbingData],
    combined: ProbingData,
    name: str = "",
) -> go.Figure:
    """Plot expected mean SNR as a function of relative total read depth."""
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
        se2 = mod_err2[None] / mod_scales[:, None, None]
        if nomod_err2 is not None:
            se2 = se2 + nomod_err2[None] / nomod_scales[:, None, None]
        se = np.sqrt(se2)
        snr = np.where(se > 0, reactivity[None] / se, 0.0)
        return np.nanmean(snr, axis=-1).mean(axis=-1)

    def _snr_band(snr_curve: np.ndarray) -> np.ndarray:
        ratio = np.where(snr_at_1 > 0, snr_curve / snr_at_1, 0.0)
        return np.abs(ratio) * snr_sem_at_1

    def _trim_leading_zeros(snr: np.ndarray) -> slice:
        nz = np.nonzero(snr > 0.01)[0]
        start = max(nz[0] - 1, 0) if len(nz) > 0 else 0
        return slice(start, None)

    xi = np.geomspace(0.1, 10, 10000)
    fig = go.Figure()

    if nomod is not None:
        mod_scales = np.maximum((xi * total_reads - nomod_reads) / mod_reads, 1e-10)
        snr_mod = _mean_snr_vec(mod_scales, np.ones_like(mod_scales))
        sem_mod = _snr_band(snr_mod)
        s = _trim_leading_zeros(snr_mod)
        xi_s = xi[s]

        fig.add_trace(
            go.Scatter(
                x=xi_s,
                y=(snr_mod - sem_mod)[s],
                line={"width": 0},
                showlegend=False,
                hoverinfo="skip",
            )
        )
        fig.add_trace(
            go.Scatter(
                x=xi_s,
                y=(snr_mod + sem_mod)[s],
                fill="tonexty",
                fillcolor=_rgba(_RDPU_070, 0.2),
                line={"width": 0},
                showlegend=False,
                hoverinfo="skip",
            )
        )
        fig.add_trace(
            go.Scatter(
                x=xi_s,
                y=snr_mod[s],
                line={"color": _RDPU_070, "width": 2},
                name="Modified",
            )
        )

        nomod_scales = np.maximum((xi * total_reads - mod_reads) / nomod_reads, 1e-10)
        snr_nomod = _mean_snr_vec(np.ones_like(nomod_scales), nomod_scales)
        sem_nomod = _snr_band(snr_nomod)
        s = _trim_leading_zeros(snr_nomod)
        xi_s = xi[s]

        fig.add_trace(
            go.Scatter(
                x=xi_s,
                y=(snr_nomod - sem_nomod)[s],
                line={"width": 0},
                showlegend=False,
                hoverinfo="skip",
            )
        )
        fig.add_trace(
            go.Scatter(
                x=xi_s,
                y=(snr_nomod + sem_nomod)[s],
                fill="tonexty",
                fillcolor=_rgba(_PUBU_070, 0.2),
                line={"width": 0},
                showlegend=False,
                hoverinfo="skip",
            )
        )
        fig.add_trace(
            go.Scatter(
                x=xi_s,
                y=snr_nomod[s],
                line={"color": _PUBU_070, "width": 2},
                name="Unmodified",
            )
        )

        curve_max = max(
            float(np.nanmax(snr_mod + sem_mod)), float(np.nanmax(snr_nomod + sem_nomod))
        )

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
        y_top = curve_max * 1.1

        fig.add_trace(
            go.Scatter(
                x=xi,
                y=pareto_upper,
                line={"width": 0},
                showlegend=False,
                hoverinfo="skip",
            )
        )
        fig.add_trace(
            go.Scatter(
                x=xi,
                y=np.full_like(xi, y_top),
                fill="tonexty",
                fillcolor="rgba(200, 200, 200, 0.3)",
                line={"width": 0},
                showlegend=False,
                hoverinfo="skip",
            )
        )
        fig.add_trace(
            go.Scatter(
                x=xi,
                y=pareto_upper,
                line={"color": "black", "width": 1, "dash": "dash"},
                name="Pareto",
            )
        )
        fig.update_yaxes(range=[0, y_top])
    else:
        snr_mod = _mean_snr_vec(xi, np.ones_like(xi))
        sem_mod = _snr_band(snr_mod)

        fig.add_trace(
            go.Scatter(
                x=xi,
                y=(snr_mod - sem_mod),
                line={"width": 0},
                showlegend=False,
                hoverinfo="skip",
            )
        )
        fig.add_trace(
            go.Scatter(
                x=xi,
                y=(snr_mod + sem_mod),
                fill="tonexty",
                fillcolor=_rgba(_RDPU_070, 0.2),
                line={"width": 0},
                showlegend=False,
                hoverinfo="skip",
            )
        )
        fig.add_trace(
            go.Scatter(
                x=xi,
                y=snr_mod,
                line={"color": _RDPU_070, "width": 2},
                name="Modified",
            )
        )

    fig.add_vline(x=1.0, line={"color": "grey", "dash": "dash", "width": 1}, opacity=0.7)

    title = f"SNR vs Read Depth ({name})" if name else "SNR vs Read Depth"
    fig.update_layout(**_base_layout(title, "Relative Total Read Depth", "Mean SNR"))
    fig.update_xaxes(type="log", range=[np.log10(0.1), np.log10(10)])
    return fig
