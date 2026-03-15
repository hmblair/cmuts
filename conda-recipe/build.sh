#!/bin/bash
set -exo pipefail

# === Setup ===

# Ensure CMake finds libraries in the conda host prefix
export CMAKE_PREFIX_PATH="${PREFIX}"

# Write version file (setuptools-scm fallback since no .git in tarball)
echo "${PKG_VERSION}" > VERSION

# Fetch htscodecs submodule (not included in GitHub release tarballs)
# NOTE: Update this commit hash when bumping the htscodecs submodule
if [[ ! -f htscodecs/htscodecs/htscodecs.h ]]; then
    rmdir htscodecs 2>/dev/null || true
    git clone https://github.com/samtools/htscodecs htscodecs
    git -C htscodecs checkout 877e6051937f85c6e5f97b70d9b6c8ab887ce81e
fi

# === Build ===

# Clean stale build artifacts (conda-build reruns for multiple Python variants,
# each with a different $PREFIX — stale CMake caches cause find_library failures)
rm -rf build
rm -rf htscodecs/build

# Build C++ binaries (htscodecs + cmake)
export HDF5_DIR="${PREFIX}"
./configure --build-only

# Verify binaries were produced
for bin in bin/cmuts-core bin/_cmuts-generate-tests; do
    if [[ ! -f "$bin" ]]; then
        echo "ERROR: $bin not produced by build"
        exit 1
    fi
done

# === Install ===

# Copy htscodecs lib into conda prefix
cp -a htscodecs/lib/libhtscodecs* "${PREFIX}/lib/"

# Install binaries
install -d "${PREFIX}/bin"
install -m 755 bin/cmuts-core "${PREFIX}/bin/"
install -m 755 bin/_cmuts-generate-tests "${PREFIX}/bin/"

# Fix htscodecs install_name and references in binaries.
# htscodecs is built with -install_name $SRC_DIR/htscodecs/.../libhtscodecs.2.dylib
# which is a build-time path that won't exist at runtime.
_INT=${INSTALL_NAME_TOOL:-install_name_tool}
if [[ "$(uname)" == "Darwin" ]]; then
    # Fix the dylib's own id
    $_INT -id "@rpath/libhtscodecs.2.dylib" "${PREFIX}/lib/libhtscodecs.2.dylib"
    # Fix references in binaries
    old_path=$(otool -L bin/cmuts-core | grep libhtscodecs | awk '{print $1}')
    if [[ -n "$old_path" ]]; then
        for b in "${PREFIX}/bin/cmuts-core" "${PREFIX}/bin/_cmuts-generate-tests"; do
            $_INT -change "$old_path" "@rpath/libhtscodecs.2.dylib" "$b"
            $_INT -add_rpath "${PREFIX}/lib" "$b" 2>/dev/null || true
        done
    fi
else
    patchelf --set-rpath "${PREFIX}/lib" "${PREFIX}/bin/cmuts-core"
    patchelf --set-rpath "${PREFIX}/lib" "${PREFIX}/bin/_cmuts-generate-tests"
fi

# Install shell scripts
for script in src/scripts/cmuts src/scripts/cmuts-align src/scripts/cmuts-generate; do
    if [[ -f "$script" ]]; then
        install -m 755 "$script" "${PREFIX}/bin/"
    fi
done

# Install Python package
export SETUPTOOLS_SCM_PRETEND_VERSION="${PKG_VERSION}"
${PYTHON} -m pip install . --no-deps --no-build-isolation -vv
