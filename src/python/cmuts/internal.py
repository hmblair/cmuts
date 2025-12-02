from __future__ import annotations
import numpy as np
import dask.array as da
import h5py
from dataclasses import dataclass
from enum import Enum
from scipy.stats import t as student
from typing import Union
import os


# ANSI escape codes (for printing)


BOLD = '\033[1m'
RESET = '\033[0m'


# Core datatypes and dataclasses


ProbingDatum = Union[np.ndarray, da.Array]

class NormScheme(Enum):
    RAW = 0
    UBR = 1
    OUTLIER = 2

class DataGroups:

    def __init__(
        self: DataGroups,
        groups: Union[list[str], None],
    ) -> None:

        if groups is None:
            groups = []

        oneD ="/counts-1d"
        twoD = "/counts-2d"

        self.oneD = [group + oneD for group in groups]
        self.twoD = [group + twoD for group in groups]

@dataclass
class Opts:
    mod: DataGroups
    nomod: DataGroups
    cutoff: int
    ins: bool
    dels: bool
    raw: bool
    outlier: bool
    blank: tuple[int, int]
    clip: tuple[bool, bool]
    sig: float


# Input array indices


IX_DEL = 4
IX_INS = 5
IX_TERM = 6


# Datasets in the output HDF5 file

@dataclass
class Datasets:
    REACTIVITY = "reactivity"
    READS = "reads"
    ERROR = "error"
    SNR = "SNR"
    NORM = "norm"
    ROI = "roi-mask"
    HEATMAP = "heatmap"
    MI = "mutual-information"
    COV = "covariance"
    SEQUENCE = 'sequence'



# Formatting


SPACE = ' '
DIV = '─'

def div() -> str:
    return SPACE * 6 + DIV * 35

def space(n: int) -> str:
    return SPACE * n

def format_s(seq: str, offset: int = 8) -> str:
    return space(offset) + seq + "\n"

def format_i(seq: str, value: float, offset: int = 8, width: int = 40) -> str:
    gap = width - len(f"{int(value):,}") - len(seq + ":") - offset
    return space(offset) + seq + ":" + space(gap) + f"{int(value):,}" + "\n"

def format_f(seq: str, value: float, offset: int = 8, width: int = 40, prec: int = 2) -> str:
    gap = width - len(f"{value:.{prec}f}") - len(seq + ":") - offset
    return space(offset) + seq + ":" + space(gap) + f"{value:.{prec}f}" + "\n"

def title(name: str, version: str) -> str:
    return space(8) + f"{name} version {version}\n" + div()


# Misc utils


def _has_arrays(
    file: h5py.File,
    groups: list[str],
) -> bool:

    if not groups:
        return False

    for group in groups:
        if group not in file:
            return False

    return True


def _accumulate_arrays(
    file: h5py.File,
    groups: list[str],
) -> da.Array:

    arr = 0

    for group in groups:
        try:
            arr += da.from_array(file[group], chunks='auto')
        except Exception as e:
            print(f"Error loading the dataset {group}:")
            raise e

    return arr


#
# Traditional MaP-seq
#


def _reactivity(
    modifications: da.Array,
    coverage: da.Array,
    mask: da.Array,
) -> da.Array:

    reactivity = da.divide(
        modifications,
        coverage,
        where=mask,
        out=da.ones_like(modifications) * np.nan,
    )
    return da.clip(reactivity, 0, 1)


def _error(
    reactivity: da.Array,
    coverage: da.Array,
    mask: da.Array,
) -> da.Array:

    return da.divide(
        da.sqrt(reactivity * (1 - reactivity)),
        da.sqrt(coverage),
        where=mask,
        out=da.ones_like(reactivity) * np.nan,
    )


def _snr(
    reactivity: da.Array,
    error: da.Array,
) -> da.Array:

    return da.divide(
        da.abs(reactivity),
        error,
        where=(error > 0),
        out=np.zeros_like(reactivity),
    ).mean(-1)


def _roi(
    seqs: int,
    seqlen: int,
    opts: Opts,
    lengths: da.Array,
) -> da.Array:

    low = opts.blank[0]
    high = opts.blank[1]

    roi = da.ones((seqs, seqlen), dtype=bool)
    if low:
        roi[:, :low] = False
    if high:
        roi[:, -high:] = False

    length_mask = (da.arange(seqlen) < lengths[:, None])
    return roi &  length_mask


#
# Pairwise MaP-seq
#


def _validate_pairwise_shape(data: da.Array) -> bool:

    return (
        data.ndim == 5 and
        data.shape[1] == data.shape[2] and
        data.shape[3] == 2 and
        data.shape[4] == 2
    )


def _probability(pairs: da.Array) -> da.Array:

    if not _validate_pairwise_shape(pairs):
        raise ValueError("The 2D counts have an incorrect shape.")

    prob = da.zeros_like(pairs)
    prob[..., 0, 0] = 1

    counts = pairs.sum((-1, -2))[..., None, None]

    return da.divide(
        pairs,
        counts,
        where=(counts > 0),
        out=prob,
    )


def _conditional(prob: da.Array) -> da.Array:

    joint = prob[..., 1, 1]
    marginal = prob[..., 1, 1] + prob[..., 1, 0]

    return da.divide(
        joint,
        marginal,
        where=(marginal > 0),
        out=np.zeros_like(joint)
    )


def _diagonal_normalize(
    values: da.Array,
    eps: float = 1E-8,
) -> da.Array:

    diag = da.diagonal(values, axis1=-2, axis2=-1)
    quot = da.sqrt(diag[..., :, None] * diag[..., None, :] + eps)
    out = da.divide(
        values,
        quot,
        where=(quot != 0),
        out=da.zeros_like(values),
    )
    return da.clip(out, -1, 1)



def _covariance(prob: da.Array) -> da.Array:

    ex = prob[..., 1, 1] + prob[..., 1, 0]
    ey = prob[..., 1, 1] + prob[..., 0, 1]
    exy = prob[..., 1, 1]

    return exy - ex * ey


def _sig_test(
    corr: da.Array,
    reads: da.Array,
    sig: float,
    eps: float = 1E-8,
) -> da.Array:

    # Bias correction

    reads = da.maximum(reads - 2, 0)

    # Calculate t-statistic and p-values

    inter = reads / (1 - corr ** 2 + eps)
    t = corr * da.sqrt(inter)
    p = 2 * (1 - da.from_array(student.cdf(da.abs(t), reads)))

    # Bonferonni correction

    n = corr.shape[-1]
    tests = n * (n - 1) // 2
    p = da.minimum(p * tests, 1.0)

    return p <= sig


def _correlation(
    prob: da.Array,
    reads: da.Array,
    p: Union[float, None],
) -> da.Array:

    cov = _covariance(prob)
    corr = _diagonal_normalize(cov)

    if p is not None:
        ix = _sig_test(corr, reads, p)
        corr = da.where(ix, corr, np.nan)

    return corr


def _mutual_information(prob: da.Array, norm: bool = True) -> da.Array:

    p_00 = prob[..., 0, 0]
    p_01 = prob[..., 0, 1]
    p_10 = prob[..., 1, 0]
    p_11 = prob[..., 1, 1]

    p_x0 = p_00 + p_01
    p_x1 = p_10 + p_11
    p_y0 = p_00 + p_10
    p_y1 = p_01 + p_11

    # p(0,0) * log(p(0,0) / (p(0) * p(0)))
    valid_00 = (p_00 > 0) & (p_x0 > 0) & (p_y0 > 0)
    ratio_00 = da.divide(
        p_00, p_x0 * p_y0,
        where=valid_00,
        out=np.ones_like(p_00)
    )
    term_00 = da.where(valid_00, p_00 * da.log(ratio_00), 0)

    # p(0,1) * log(p(0,1) / (p(0) * p(1)))
    valid_01 = (p_01 > 0) & (p_x0 > 0) & (p_y1 > 0)
    ratio_01 = da.divide(
        p_01, p_x0 * p_y1,
        where=valid_01,
        out=np.ones_like(p_01)
    )
    term_01 = da.where(valid_01, p_01 * da.log(ratio_01), 0)

    # p(1,0) * log(p(1,0) / (p(1) * p(0)))
    valid_10 = (p_10 > 0) & (p_x1 > 0) & (p_y0 > 0)
    ratio_10 = da.divide(
        p_10, p_x1 * p_y0,
        where=valid_10,
        out=np.ones_like(p_10)
    )
    term_10 = da.where(valid_10, p_10 * da.log(ratio_10), 0)

    # p(1,1) * log(p(1,1) / (p(1) * p(1)))
    valid_11 = (p_11 > 0) & (p_x1 > 0) & (p_y1 > 0)
    ratio_11 = da.divide(
        p_11, p_x1 * p_y1,
        where=valid_11,
        out=np.ones_like(p_11)
    )
    term_11 = da.where(valid_11, p_11 * da.log(ratio_11), 0)

    mi = term_00 + term_01 + term_10 + term_11
    if norm:
        return _diagonal_normalize(mi)
    else:
        return mi


#
# Main normalization function
#


def _pairwise_mask(mask: da.Array) -> da.Array:

    return mask[..., :, None] & mask[..., None, :]


def _data_from_counts(
    counts: da.Array,
    pairs: Union[da.Array, None],
    opts: Opts,
    lengths: da.Array
) -> tuple[Union[da.Array, None], ...]:

    seqs = int(counts.shape[0])
    seqlen = int(counts.shape[1])

    roi = _roi(seqs, seqlen, opts, lengths)

    # Compute the modification heatmap
    # The termination events are normalized separately

    heatmap = da.mean(counts, axis=(0, 1))
    heatmap = heatmap / heatmap.sum()

    # Compute the coverage as the sum of all matches and mismatches.
    # Store the mean coverage as a function of position for later.

    coverage = da.sum(counts[..., :IX_DEL], axis=(2, 3))
    mc = da.mean(coverage, axis=0)

    # Compute the coverage mask

    mask = (coverage >= opts.cutoff)
    mask[~roi] = False

    # Compute termination events

    term = da.sum(counts[..., IX_TERM], axis=2)

    # Compute the trace to get the non-modified rates only

    trace = da.trace(counts, axis1=2, axis2=3)

    # Subtract the trace to get the pure modifications

    modifications = coverage - trace

    # Account for insertions and deletions

    if opts.ins:
        modifications += da.sum(counts[..., IX_INS], axis=-1)
    if opts.dels:
        modifications += da.sum(counts[..., IX_DEL], axis=-1)

    # Divide by the coverage to get the un-normalised reactivity

    reactivity = _reactivity(modifications, coverage, mask)

    # SEM of a Bernoulli random variable

    err = _error(reactivity, coverage, mask)

    # Get the read counts and compute all values

    reads = da.max(coverage, axis=(-1,))

    # Signal-to-noise

    snr = _snr(reactivity, err)

    # Pairwise probabilities

    if pairs is not None:

        prob = _probability(pairs)

        mask_2d = _pairwise_mask(mask)

        pairs = pairs.sum((-2, -1))
        covariance = _correlation(prob, pairs, opts.sig)
        mi = _mutual_information(prob)

        covariance[~mask_2d] = np.nan
        mi[~mask_2d] = np.nan

    else:
        prob = None
        covariance = None
        mi = None

    return (
        reactivity,
        reads,
        err,
        snr,
        mask,
        heatmap,
        mc,
        term,
        pairs,
        prob,
        covariance,
        mi,
    )


def _print_stats_single(
    data: ProbingData,
) -> str:

    multi = (data.reactivity.shape[0] > 1)

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
    with np.errstate(divide='ignore', invalid='ignore'):
        med2m = mr / medr

    out = ""

    if multi:
        out += format_f("Mean reads", mr, 10)
        out += format_i(f"Median reads", medr, 10)
        out += format_f(f"Mean-to-Median", med2m, 10)
        out += format_f("Mean reactivity", mrr, 10, prec=3)
        out += format_f("Mean error", me, 10, prec=3)
        out += format_f("Deletion rate", del_r, 10, prec=3)
        out += format_f("Dropout fraction", dropout, 10)
        out += format_f("Mean SNR", msnr, 10)
        out += format_f("SNR > 1", fsnr, 10)
    else:
        out += format_i("Reads", mr, 10)
        out += format_f("Mean reactivity", mrr, 10, prec=3)
        out += format_f("Mean error", me, 10, prec=3)
        out += format_f("Deletion rate", del_r, 10, prec=3)
        out += format_f("SNR", msnr, 10)

    return out


class ProbingData:

    def __init__(
        self: ProbingData,
        sequences: Union[da.Array, None],
        reactivity: ProbingDatum,
        reads: ProbingDatum,
        error: ProbingDatum,
        snr: ProbingDatum,
        mask: ProbingDatum,
        heatmap: ProbingDatum,
        coverage: ProbingDatum,
        terminations: ProbingDatum,
        pairs: Union[ProbingDatum, None],
        probability: Union[ProbingDatum, None],
        covariance: Union[ProbingDatum, None],
        mi: Union[ProbingDatum, None],
    ) -> None:

        self.sequences = sequences
        self.norm = None

        # 1D values

        self.reactivity = reactivity
        self.reads = reads
        self.error = error
        self.snr = snr
        self.mask = mask
        self.heatmap = heatmap
        self.coverage = coverage
        self.terminations = terminations

        # 2D values

        self.pairs = pairs
        self.probability = probability
        self.covariance = covariance
        self.mi = mi

    def size(
        self: ProbingData,
    ) -> int:

        return int(self.reads.shape[0])

    def single(
        self: ProbingData,
    ) -> bool:

        return self.size() == 1

    def __repr__(
        self: ProbingData,
    ) -> str:

        return _print_stats_single(self)

    def compute(
        self: ProbingData,
    ) -> ProbingData:

        return ProbingData(
            *da.compute(
                self.sequences,
                self.reactivity,
                self.reads,
                self.error,
                self.snr,
                self.mask,
                self.heatmap,
                self.coverage,
                self.terminations,
                self.pairs,
                self.probability,
                self.covariance,
                self.mi,
            )
        )

    def clip(
        self: ProbingData,
        low: bool,
        high: bool,
    ) -> None:

        if low:
            self.reactivity = da.clip(self.reactivity, 0, None)
        if high:
            self.reactivity = da.clip(self.reactivity, None, 1)
        if low or high:
            self.snr = _snr(self.reactivity, self.error)

    def normalize(
        self: ProbingData,
        norm: ProbingDatum,
    ) -> None:

        norm = da.nan_to_num(norm, 1)

        if norm.ndim > 0:
            norm[norm <= 0] = 1
        elif norm.item() <= 0:
            norm = 1

        self.reactivity /= norm
        self.error /= norm
        self.norm = norm

    def save(
        self: ProbingData,
        group: str,
        file: str,
    ) -> None:

        data = {
            group + '/' + Datasets.REACTIVITY: da.from_array(self.reactivity),
            group + '/' + Datasets.READS: da.from_array(self.reads),
            group + '/' + Datasets.NORM: da.from_array(self.norm),
            group + '/' + Datasets.ERROR: da.from_array(self.error),
            group + '/' + Datasets.SNR: da.from_array(self.snr),
            group + '/' + Datasets.HEATMAP: da.from_array(self.heatmap),
        }

        if self.mi is not None:
            data[group + '/' + Datasets.MI] = da.from_array(self.mi)

        if self.covariance is not None:
            data[group + '/' + Datasets.COV] = da.from_array(self.covariance)

        if self.sequences is not None:
            data[Datasets.SEQUENCE] = da.from_array(self.sequences)

        da.to_hdf5(file, data)

    @classmethod
    def load(
        cls,
        group: str,
        file: h5py.File,
    ) -> 'ProbingData':
        """Load ProbingData from an HDF5 file."""

        def get(name: str) -> np.ndarray:
            key = f"{group}/{name}" if group else name
            return np.array(file[key]) if key in file else None

        sequences = np.array(file[Datasets.SEQUENCE]) if Datasets.SEQUENCE in file else None
        reactivity = get(Datasets.REACTIVITY)
        reads = get(Datasets.READS)
        error = get(Datasets.ERROR)
        snr = get(Datasets.SNR)
        heatmap = get(Datasets.HEATMAP)
        mi = get(Datasets.MI)
        covariance = get(Datasets.COV)
        norm = get(Datasets.NORM)

        # Create minimal ProbingData with required fields
        # Fields not saved are set to empty arrays or None
        n = reactivity.shape[0] if reactivity is not None else 0
        length = reactivity.shape[1] if reactivity is not None and reactivity.ndim > 1 else 0

        data = cls(
            sequences=sequences,
            reactivity=reactivity,
            reads=reads,
            error=error,
            snr=snr,
            mask=np.ones_like(reactivity, dtype=bool) if reactivity is not None else None,
            heatmap=heatmap,
            coverage=np.sum(~np.isnan(reactivity), axis=0) if reactivity is not None else None,
            terminations=np.zeros((n, length)) if reactivity is not None else None,
            pairs=None,
            probability=None,
            covariance=covariance,
            mi=mi,
        )
        data.norm = norm
        return data


def _sequences_from_counts(
    file: h5py.File,
    dataset: str = "sequences",
) -> Union[da.Array, None]:

    # Get the tokenized sequences if they exist

    if dataset in file:
        return da.from_array(file[dataset], chunks='auto')

    return None


def _probing_data_from_counts(
    file: h5py.File,
    groups: DataGroups,
    opts: Opts,
    lengths: da.Array
    ) -> Union[ProbingData, None]:

        if not groups.oneD:
            return None

        sequences = _sequences_from_counts(file)

        counts = _accumulate_arrays(file, groups.oneD)
        if _has_arrays(file, groups.twoD):
            pairs = _accumulate_arrays(file, groups.twoD)
        else:
            pairs = None

        data = _data_from_counts(counts, pairs, opts, lengths)
        return ProbingData(sequences, *data)


def _merge_conditions(
    mod: ProbingData,
    nomod: Union[ProbingData, None],
) -> ProbingData:

    if nomod is None:
        return mod

    sequences = mod.sequences

    # 1D data

    reactivity = mod.reactivity - nomod.reactivity
    reads = mod.reads + nomod.reads
    error = np.sqrt(mod.error ** 2 + nomod.error ** 2)
    snr = _snr(reactivity, error)
    mask = mod.mask & nomod.mask
    heatmap = mod.heatmap
    coverage = mod.coverage + nomod.coverage
    terminations = mod.terminations + nomod.terminations

    # 2D data

    pairs = mod.pairs
    probability = mod.probability
    mi = mod.mi
    covariance = mod.covariance

    return ProbingData(
        sequences,
        reactivity,
        reads,
        error,
        snr,
        mask,
        heatmap,
        coverage,
        terminations,
        pairs,
        probability,
        covariance,
        mi,
    )


def _get_norm_scheme(raw: bool, outlier: bool) -> NormScheme:

    if raw and outlier:
        raise ValueError("The --raw and --norm-outlier flags are mutually exclusive.")

    if raw:
        return NormScheme.RAW
    if outlier:
        return NormScheme.OUTLIER
    else:
        return NormScheme.UBR


def _get_norm_raw(
    data: ProbingData,
    opts: Opts,
) -> np.ndarray:
    """
    No normalization.
    """

    return np.ones(1, dtype=data.reactivity.dtype)


def _get_norm_percentile(
    data: ProbingData,
    opts: Opts,
) -> np.ndarray:

    READ_CUTOFF = 100
    PERCENTILE = 90

    _good_pos = data.mask & (data.reads > READ_CUTOFF)[:, None]
    high_reactivity = data.reactivity[_good_pos]

    if not high_reactivity.size:
        return _get_norm_raw(data, opts)

    return np.percentile(high_reactivity.flatten(), PERCENTILE)


def _get_norm_outlier(
    data: ProbingData,
    opts: Opts,
) -> np.ndarray:
    """
    2-8% Normalization: From top 10%, ignore top 2% and divide by average of
    remaining 8%.
    """

    def norm_1d(
        arr: np.ndarray,
        low: float = 0.02,
        high: float = 0.1,
    ) -> float:

        arr = arr[~np.isnan(arr)]

        p2 = max(1, round(len(arr) * low) - 1)
        p10 = max(1, round(len(arr) * high) - 1)

        sarr = np.sort(arr)[::-1]
        if sarr.size == 0:
            m = 1
        else:
            m = np.mean(sarr[p2:p10+1])

        return float(m)

    return np.apply_along_axis(norm_1d, axis=1, arr=data.reactivity)


def _get_norm(
    data: ProbingData,
    opts: Opts,
) -> np.ndarray:

    scheme = _get_norm_scheme(opts.raw, opts.outlier)

    if scheme == NormScheme.RAW:
        return _get_norm_raw(data, opts)
    if scheme == NormScheme.UBR:
        return _get_norm_percentile(data, opts)
    if scheme == NormScheme.OUTLIER:
        return _get_norm_outlier(data, opts)

    raise ValueError(f"Invalid normalization scheme {scheme}.")


def _lengths_from_fasta(
    file: str
) -> da.Array:

    lengths = []
    length = 0
    start = True

    with open(file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                if start:
                    start = False
                    continue
                lengths.append(length)
                length = 0
            else:
                length += len(line.strip())

    lengths.append(length)
    return da.array(lengths)


def compute_reactivity(
    file: h5py.File,
    fasta: str,
    opts: Opts,
) -> tuple[ProbingData, Union[ProbingData, None], ProbingData]:
    """
    Compute reactivity profiles from mutation count data.

    Args:
        file: Open HDF5 file containing count data
        fasta: Path to FASTA file with reference sequences
        opts: Processing options

    Returns:
        Tuple of (modified, unmodified, combined) ProbingData objects

    Raises:
        FileNotFoundError: If FASTA file does not exist
    """
    if not os.path.exists(fasta):
        raise FileNotFoundError(f"FASTA file not found: {fasta}")

    # Get reactivity and associated values

    lengths = _lengths_from_fasta(fasta)
    mod = _probing_data_from_counts(file, opts.mod, opts, lengths)
    nomod = _probing_data_from_counts(file, opts.nomod, opts, lengths)

    assert mod is not None
    combined = _merge_conditions(mod, nomod)

    # Compute

    mod = mod.compute()
    combined = combined.compute()
    if nomod is not None:
        nomod = nomod.compute()

    # Normalize

    norm = _get_norm(combined, opts)

    mod.normalize(norm)
    combined.normalize(norm)
    if nomod is not None:
        nomod.normalize(norm)

    # Clip

    mod.clip(opts.clip[0], opts.clip[1])
    combined.clip(opts.clip[0], opts.clip[1])
    if nomod is not None:
        nomod.clip(opts.clip[0], opts.clip[1])

    return mod, nomod, combined


def stats(
    mod: ProbingData,
    nomod: Union[ProbingData, None],
    combined: ProbingData,
) -> None:

    multi = (combined.reactivity.shape[0] > 1)

    nseqs = combined.reads.shape[0]
    seqlen = combined.reactivity.shape[1]
    reads = combined.reads.sum()

    print(format_i("References", nseqs, 8), end="")
    print(format_i("Reference length", seqlen, 8), end="")
    if multi:
        print(format_i("Total reads", reads, 8), end="")

    print(format_s(f"{BOLD}Treated (MOD):{RESET}", 8), end="")
    print(mod, end="")

    if nomod is not None:
        print(format_s(f"{BOLD}Untreated (NOMOD):{RESET}", 8), end="")
        print(nomod, end="")
        print(format_s(f"{BOLD}Combined:{RESET}", 8), end="")
        print(combined, end="")

    print(div())
    print()
