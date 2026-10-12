# Environment and Builds

## Setting up

The tests need the build dependencies in the README, `samtools` on `PATH`, and a Python environment with the test dependencies. The tests of `cmuts align` also need `minimap2`, and `vsearch` for the paired-end case. Without them, those tests are skipped.

```sh
uv venv .venv
uv pip install --python .venv/bin/python --group dev
```

Rendering the documentation also needs the `docs` group:

```sh
uv pip install --python .venv/bin/python --group docs
```

## Building

```sh
make
```

Each variant builds into its own directory under `build/`. `CC` sets the compiler and `OPT` sets the optimization flags.

The build finds headers and libraries through `pkg-config` where htslib and HDF5 provide a `.pc` file, and through plain `-l` and `CPATH` otherwise.

## Tests

```sh
make check
```

This builds the programs and runs pytest over `tests/` against the binaries in the build directory. Pass `PYTHON` to `make` to use a Python other than `.venv/bin/python`.

The tests refuse a binary outside the repository. To run pytest directly, put the build directory first on `PATH`:

```sh
PATH=$PWD/build/release:$PATH .venv/bin/python -m pytest tests/test_filtering.py -x
```

## Linter

```sh
make lint
```

This writes `compile_commands.json` and runs clang-tidy over every source file under `src/` and `apps/` with the checks in `.clang-tidy`. Set `CLANG_TIDY` to use a different linter binary.

## Sanitizers

`SAN` names a sanitizer to build with, and applies to the target that follows it:

```sh
make SAN=asan         # address and undefined behavior
make check SAN=tsan   # thread, over the tests
```

The tests capture the standard error of every program they run, and a failure reports only the exit status. `make check` therefore writes each sanitizer report to its own file under the build directory of the variant:

```
build/tsan/logs/tsan.<pid>
```

A value already set in `ASAN_OPTIONS` or `TSAN_OPTIONS` overrides the one `make` sets.

## Documentation

The documentation of the program arguments and of the HDF5 output is generated from the binaries.

```sh
make docs   # rewrite the generated blocks under docs/
make site   # additionally render the site into site/
```

When the arguments or the output format of a program change, run `make docs` and include its changes in the commit.

The site holds one rendering per release from `v2.0.0` on, each under its tag, and a rendering of the working tree under `latest`. A dropdown in the sidebar switches between them. `make site` reads the docs tree of each release from its tag, so the tags must be present locally, but nothing needs to be checked out. The root of the site redirects to `latest`.

To preview the site locally:

```sh
make serve   # serve site/ at http://localhost:8000/
```

`PORT` chooses another port. The target refuses to start when the port is already in use.

## Static Binaries

`scripts/static-deps.sh` builds htslib, HDF5, and their compression libraries as static archives under one prefix. The sources are pinned by version and checksum.

```sh
scripts/static-deps.sh deps-static
make STATIC_PREFIX=deps-static
```

On macOS, the system libraries stay dynamic. A musl compiler, passed as the second argument of `scripts/static-deps.sh` as in `scripts/static-deps.sh deps-static musl-gcc`, or any Alpine gcc, makes the whole link static. `make check STATIC_PREFIX=deps-static` runs the tests against the result. `scripts/package.sh` wraps the result in the tarball that a release carries.

# Contributing

## Writing Tests

A test that runs a program over an alignment takes the `data` fixture. The fixture is parametrized over the alignments in `tests/datasets.py` and over every format htslib reads, and yields one `Dataset` per combination. A new test must take this fixture and must test a contract that holds for every input alignment and every format. The `catalogue` fixture gives a lookup by name and format, for the few tests that need a second dataset.

A test that needs specific formats pins `fmt` instead of taking all three. A module pins it by overriding the fixture, and a single test pins it with `indirect`.

The `falsifiable` fixture records whether the test can fail on the current dataset, which keeps tests from passing vacuously. Call it between the setup and the assertions of the test. A test that no dataset can falsify fails the run, which shows that `tests/datasets.py` needs a new dataset.

## Continuous Integration

Two workflows run on every push to `main` and on every pull request against it.

`CI` runs `make check` with gcc and clang on Linux and with clang on macOS, and `make lint` on Linux. It does not run the sanitizers, which take longer than the rest of the run together. Run those locally.

On a push to `main`, after those jobs pass, `CI` pins the `DasLab/cmuts-space` Dockerfile to the commit. It then waits for the space to rebuild and redeploy, and fails if that fails, so a change that breaks the space is reported here.

`Documentation` renders the site, and fails if the generated docs are not current. It does not update them. A push to `main` also deploys the rendered site.

## Releasing

A release is a `v*` tag on `main`. Before tagging, set three files to the version the tag will name:

1. `include/version.h`, which the programs report and which names the tarballs.
2. `CITATION.cff`, whose `version` and `date-released` fields a citation of the software carries. Zenodo archives each tag and reads this file for the title, the authors, and the license of the record.
3. `CHANGELOG.md`, which needs an entry under `## [<version>]`. The release takes its notes from that entry.

Pushing the tag runs the `Release` workflow. Its first job runs `scripts/check-tag.sh`, which fails if the tag disagrees with the header or the citation file, and `scripts/release-notes.sh`, which fails if the changelog holds no entry for the version. The platform jobs wait on that job, so a tag that disagrees with any of the three files fails before anything is built.

The platform jobs build the static binaries for Linux (x86_64 and aarch64, in Alpine containers) and macOS (arm64), run the test suite against each, and attach the `scripts/package.sh` tarballs to a draft GitHub release with the changelog entry as its notes. You publish the draft by hand.

The tarballs are named `cmuts-<version>-<os>-<arch>.tar.gz`, where `<os>` is `uname -s` in lower case and `<arch>` is `uname -m`. The release page lists them only after its scripts run, so a fetch of the page from a terminal shows no assets. Each is served at a fixed URL that a script can build from the pattern:

```sh
curl -LO https://github.com/DasLab/cmuts/releases/download/v<version>/cmuts-<version>-$(uname -s | tr '[:upper:]' '[:lower:]')-$(uname -m).tar.gz
```
