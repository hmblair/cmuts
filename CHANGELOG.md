# Changelog

## [1.2.4] - 2025-01-22

### Added

- Configurable logging system via `CMUTS_LOG_LEVEL` and `CMUTS_LOG_STDERR` environment variables
- `TRACE` log level for detailed MPI debugging

### Fixed

- **Critical**: MPI deadlock in HDF5 collective I/O when processes had asymmetric workloads
- Root cause: `skip` flag was reset after each read, causing processes to participate in different numbers of collective writes

### Changed

- Version now derived from git tags (configure generates VERSION file)

## [1.2.3] - 2025-12-01

### Added

- Pytest integration test suite (`cmuts test`) with edge case, CRAM, fuzzing, and stress tests
- `--verbose` flag to `cmuts core` for debug logging
- Structured logging with levels (DEBUG, INFO, WARN, ERROR) to `.cmuts.log`
- `Alignment::create()` factory method to prevent field order initialization bugs
- Validation macros (`CMUTS_CHECK_BOUNDS`, `CMUTS_CHECK_NOT_NULL`, etc.)
- Linux installation instructions (Ubuntu/Debian, Fedora/RHEL) to documentation
- `cmuts test` documentation page

### Fixed

- **Critical**: CRAM field order bug where `primary` and `reversed` fields were swapped, causing all CRAM reads to be filtered as non-primary
- Linux compatibility in profile script (portable GNU time detection)

### Changed

- `cmuts test` now runs pytest tests by default (use `--quick` for fast tests, `--cram` for CRAM-only)
- Documentation: fixed outdated `cmuts pairwise` reference in pipeline.md

### Removed

- Legacy bash-based test script (`cmuts-test`, `cmuts-test-compare`)

## [1.2.2] - 2025-11-19

### Added

- Fast-path for processing sparse alignments (reads << references)
- `--secondary` flag for processing secondary alignments
- JSON-based argument file for `cmuts core` and `cmuts generate`

### Changed

#### Non-Breaking

- Minor bugfixes
- Reads with missing mapping qualities are skipped, unless `--min-mapq` is set to `0`
- Seek BAM/CRAM only on MPI builds for slightly faster processing
- Replaced file-based mutex with semaphores

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
