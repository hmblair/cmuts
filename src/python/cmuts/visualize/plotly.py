"""Plotly implementations of cmuts plots.

Each function returns a ``go.Figure`` for interactive display (e.g. in Gradio).
The matplotlib module (``plotting.py``) remains the CLI default and saves PNGs.
"""

from typing import Any, Optional

import numpy as np
import plotly.colors as pc
import plotly.graph_objects as go

from cmuts.internal import SNRCurves

from . import _transforms

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


def _add_band(
    fig: go.Figure,
    x: np.ndarray,
    lower: np.ndarray,
    upper: np.ndarray,
    fillcolor: str,
) -> None:
    """Shade the region between ``lower`` and ``upper``: an invisible lower
    trace followed by an upper trace that fills down to it (``fill="tonexty"``).
    """
    fig.add_trace(go.Scatter(x=x, y=lower, line={"width": 0}, showlegend=False, hoverinfo="skip"))
    fig.add_trace(
        go.Scatter(
            x=x,
            y=upper,
            fill="tonexty",
            fillcolor=fillcolor,
            line={"width": 0},
            showlegend=False,
            hoverinfo="skip",
        )
    )


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

    # Intentionally not using _base_layout: heatmaps need square cells
    # (yaxis scaleanchor + constrain "domain"), a reversed y-axis, and fixed
    # sizing that the shared line-plot layout does not provide. Keep in sync
    # with _base_layout's font sizing by hand if the base style changes.
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
    lr = _transforms.read_log_depth(reads)

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


def plot_reads_per_block(reads: np.ndarray, name: str = "", nblocks: int = 100) -> go.Figure:
    """Plot mean reads per reference, with the references split into a
    fixed number of equal-sized blocks in their FASTA order.

    The x-axis is the block index (0..nblocks-1); the y-axis is the
    arithmetic mean of per-reference read counts inside that block. No
    sorting is applied, so spatial trends along the FASTA are preserved.
    """
    means = _transforms.reads_per_block(reads, nblocks)

    title = f"Mean reads per reference bin ({name})" if name else "Mean reads per reference bin"
    fig = go.Figure()
    _line_trace(fig, means)
    fig.update_layout(**_base_layout(title, "Reference bin", "Mean reads"))
    fig.update_xaxes(range=[0, len(means)])
    return fig


def plot_termination(term: np.ndarray, name: str = "") -> go.Figure:
    term = _transforms.termination_density(term)

    title = f"Termination by position ({name})" if name else "Termination by position"
    fig = go.Figure()
    _line_trace(fig, term)
    fig.update_layout(**_base_layout(title, "Residue", "Termination density"))
    fig.update_xaxes(range=[0, len(term)])
    fig.update_yaxes(range=[0, 1.05])
    return fig


def _seq_tick_labels(sequence: str) -> tuple[list[int], list[str]]:
    """Build x-axis tick values/labels showing nucleotide letters.

    Each position gets a letter; positions 1, every 10th, and the last also
    get the residue number on a second line below the letter.
    """
    n = len(sequence)
    tickvals = list(range(n))
    ticktext: list[str] = []
    for i, nt in enumerate(sequence):
        pos = i + 1
        if pos == 1 or pos % 10 == 0 or pos == n:
            ticktext.append(f"{nt}<br>{pos}")
        else:
            ticktext.append(nt)
    return tickvals, ticktext


def _apply_sequence_axis(fig: go.Figure, sequence: Optional[str], length: int) -> None:
    """Replace the x-axis numeric ticks with sequence letter ticks (if seq fits)."""
    if not sequence or len(sequence) < length:
        return
    tickvals, ticktext = _seq_tick_labels(sequence[:length])
    fig.update_xaxes(tickmode="array", tickvals=tickvals, ticktext=ticktext)


def plot_profile(
    reactivity: np.ndarray,
    error: Optional[np.ndarray] = None,
    name: str = "",
    sequence: Optional[str] = None,
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
    _apply_sequence_axis(fig, sequence, len(reactivity))
    return fig


def plot_profiles(
    reactivities: list[np.ndarray],
    names: list[str],
    sequence: Optional[str] = None,
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
    _apply_sequence_axis(fig, sequence, max_len)
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
    sequence: Optional[str] = None,
) -> go.Figure:
    reactivity = np.asarray(reactivity)
    error = np.asarray(error)
    if reactivity.shape[0] == 1:
        return plot_profile(reactivity[0], error[0], name, sequence=sequence)
    else:
        return plot_multiple_examples(reactivity[:num], name)


def plot_coverage(coverage: np.ndarray, reads: np.ndarray, name: str = "") -> go.Figure:
    data = _transforms.coverage_fraction(coverage, reads)

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
    vlow = _transforms.pairwise_coverage_bounds(values).vlow
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
    fig.update_layout(yaxis_autorange="reversed", yaxis_scaleanchor="x")
    return fig


def plot_correlation(values: np.ndarray, name: str = "") -> go.Figure:
    linthresh = _transforms.CORRELATION_LINTHRESH
    transformed = _transforms.symlog(values, linthresh)
    t_max = _transforms.symlog(np.array([1.0]), linthresh)[0]

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
                    _transforms.symlog(np.array([v]), linthresh)[0]
                    for v in [-1, -0.1, -0.01, 0, 0.01, 0.1, 1]
                ],
                "ticktext": ["-1", "-0.1", "-0.01", "0", "0.01", "0.1", "1"],
            },
        )
    )
    fig.update_layout(**_base_layout(title, "Residue", "Residue"))
    fig.update_layout(yaxis_autorange="reversed", yaxis_scaleanchor="x")
    return fig


def plot_mi(values: np.ndarray, name: str = "") -> go.Figure:
    vlow, vhigh = _transforms.mi_bounds(values)

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
    fig.update_layout(yaxis_autorange="reversed", yaxis_scaleanchor="x")
    return fig


def plot_snr_scaling(curves: SNRCurves, name: str = "") -> go.Figure:
    """Render precomputed SNR-vs-read-depth curves (see cmuts.compute_snr_curves)."""

    def _trim(snr: np.ndarray) -> slice:
        nz = np.nonzero(snr > 0.01)[0]
        start = max(nz[0] - 1, 0) if len(nz) > 0 else 0
        return slice(start, None)

    xi = curves.xi
    fig = go.Figure()

    if curves.nomod is not None:
        # nomod_sem, pareto and pareto_sem are populated together with nomod.
        assert curves.nomod_sem is not None
        assert curves.pareto is not None and curves.pareto_sem is not None
        s = _trim(curves.mod)
        _add_band(
            fig,
            xi[s],
            (curves.mod - curves.mod_sem)[s],
            (curves.mod + curves.mod_sem)[s],
            _rgba(_RDPU_070, 0.2),
        )
        fig.add_trace(
            go.Scatter(
                x=xi[s], y=curves.mod[s], line={"color": _RDPU_070, "width": 2}, name="Modified"
            )
        )

        s = _trim(curves.nomod)
        _add_band(
            fig,
            xi[s],
            (curves.nomod - curves.nomod_sem)[s],
            (curves.nomod + curves.nomod_sem)[s],
            _rgba(_PUBU_070, 0.2),
        )
        fig.add_trace(
            go.Scatter(
                x=xi[s], y=curves.nomod[s], line={"color": _PUBU_070, "width": 2}, name="Unmodified"
            )
        )

        curve_max = max(
            float(np.nanmax(curves.mod + curves.mod_sem)),
            float(np.nanmax(curves.nomod + curves.nomod_sem)),
        )
        y_top = curve_max * 1.1
        pareto_upper = curves.pareto + curves.pareto_sem
        _add_band(fig, xi, pareto_upper, np.full_like(xi, y_top), "rgba(200, 200, 200, 0.3)")
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
        _add_band(
            fig, xi, curves.mod - curves.mod_sem, curves.mod + curves.mod_sem, _rgba(_RDPU_070, 0.2)
        )
        fig.add_trace(
            go.Scatter(x=xi, y=curves.mod, line={"color": _RDPU_070, "width": 2}, name="Modified")
        )

    fig.add_vline(x=1.0, line={"color": "grey", "dash": "dash", "width": 1}, opacity=0.7)
    title = f"SNR vs Read Depth ({name})" if name else "SNR vs Read Depth"
    fig.update_layout(**_base_layout(title, "Relative Total Read Depth", "Mean SNR"))
    fig.update_xaxes(type="log", range=[np.log10(xi[0]), np.log10(xi[-1])])
    return fig
