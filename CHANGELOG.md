# Changelog

## [1.2.1] - 2025-11-14

### Added

- New subprogram `cmuts plot` to plot multiple reactivity profiles, or profiles generated outside of `cmuts normalize`
- Re-instated `--print-every` argument of `cmuts core`, now measured in seconds rather than reads
- Added `--downsample` argument to `cmuts core`, to limit read depths per reference
- Lock for when multiple processes attempt to generate `.cmfa` or `.cmix` files

### Changed

#### Non-Breaking

- Upgraded algorithm for detecting and spreading ambiguous deletions
- `cmuts normalize` now compatible with Python 3.9 and Python 3.16

## [1.2.0] - 2025-11-07

### Added

- Strand filtering via `--no-reverse` and `--only-reverse`
- Hamming distance filtering via `--max-hamming`
- Base filtering via `--ignore-bases`, such as for DMS probing
- Speedup of pairwise profile processing via OMP
- `Low-quality bases` stat added to `cmuts core`
- Length-based blanking for libraries with variable-length references in `cmuts normalize`
- Additional error and Bonferroni-corrected correlation plots in `cmuts normalize`
- Optional `--sig` flag in `cmuts normalize` to specify correlation significance levels
- `--no-insertions` and `--no-deletions` added to `cmuts normalize`
- Python API to `cmuts visualize` (`cmuts.visualize`) and `cmuts normalize` (`cmuts.normalize`)
- `--clean` flag for `configure` script to build from scratch

### Changed

#### Breaking

- `--joint` argument of `cmuts core` renamed `--pairwise`
- Output datasets from `cmuts core` now named `counts-1d` and `counts-2d` for normal and `--pairwise` mode, respectively
- Merged `cmuts pairwise` into `cmuts normalize`, which will auto-process 2D counts if present alongside 1D counts
- `cmuts normalize` now requires a `--fasta` input

#### Non-Breaking

- 1D counts computed even when running `cmuts core` in `--pairwise` mode
- `configure` script will `pip`-install the `cmuts` Python API

### Removed

- Unsused `--subsample` functionality
- `--print-every` argument (now time-based printing)
- Unsused and only partially supported `--low-mem` runtime mode
- Now-merged `cmuts pairwise`
