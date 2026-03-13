# cmuts

Chemical mutation counting for RNA MaP-seq experiments. C++ core with OpenMP, HDF5 output, wrapped by a `cmuts` bash dispatcher.

## Build

```bash
./configure            # build C++ binaries, generate docs, pip install
./configure --clean    # full clean rebuild
./configure --mpi      # build with MPI support
./configure --debug    # build with sanitizers
./configure --install-completions  # install bash/zsh completions
```

Requires: cmake, make, libhdf5, libhts, zlib. On macOS: `brew install hdf5 htslib zlib cmake`.

## Test

```bash
cmuts test             # run full test suite
cmuts test --quick     # skip slow tests
cmuts test --cram      # CRAM-specific tests only
```

## Project structure

- `src/cpp/` — C++ source (core mutation counter, test data generator)
- `src/python/` — Python CLI scripts (normalize, plot, visualize)
- `src/scripts/` — Bash CLI scripts (cmuts dispatcher, cmuts-align)
- `config/` — JSON argument definitions (used to generate C++ headers, docs, and completions)
- `scripts/` — Code generators (args, docs, completions)
- `tests/` — Python integration tests
