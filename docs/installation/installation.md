# Basic Installation

This page covers the installation of the single-threaded version of `cmuts`. Before starting, ensure you have all [required dependencies](requirements.md) installed on your system.

!!! tip "Verify Dependencies"
    Run these commands to verify your system is ready:
    ```bash
    cmake --version               # Should be 3.29+
    samtools --version            # Should show version info
    autoreconf --version          # Should show version info
    pkg-config --modversion hdf5  # Should show HDF5 version
    ```

## Installation

### 1. Obtain `cmuts`

Clone the `cmuts` repository with all submodules:

```bash
git clone --recurse-submodules https://github.com/hmblair/cmuts
cd cmuts
```


### 2. Build `cmuts`

The configuration script contains all required setup and build steps.

```
./configure
```

If building the multithreaded version, pass the `--mpi` flag to the script.

```
./configure --mpi
```

### 3. Modify PATH

Afterwards, add the `cmuts` binary directory to your PATH:

=== "Temporary (current session)"
    ```bash
    export PATH="$(pwd)/bin:$PATH"
    ```

=== "Permanent (bash)"
    ```bash
    echo 'export PATH="'$(pwd)'/bin:$PATH"' >> ~/.bashrc
    source ~/.bashrc
    ```

=== "Permanent (zsh)"
    ```bash
    echo 'export PATH="'$(pwd)'/bin:$PATH"' >> ~/.zshrc
    source ~/.zshrc
    ```

## Troubleshooting

### Common Build Issues

!!! failure "cmake cannot find HDF5"
    **Solution**: Set the HDF5_DIR environment variable
    ```bash
    # For Homebrew on macOS
    export HDF5_DIR=$(brew --prefix hdf5)
    ./configure
    
    # For Ubuntu/Debian
    export HDF5_DIR=/usr/lib/x86_64-linux-gnu/hdf5/serial
    ./configure
    ```

!!! failure "samtools not found"
    **Solution**: Ensure samtools is in your PATH
    ```bash
    which samtools  # Should return a path
    # If not found, install via package manager or add to PATH
    ```

!!! failure "autoconf missing"
    **Solution**: Install autoconf
    ```bash
    # macOS
    brew install autoconf

    # Ubuntu/Debian
    sudo apt install autoconf
    ```

!!! failure "submodules not initialized"
    **Error**: `fatal: not a git repository` during build

    **Solution**: Ensure you cloned with submodules
    ```bash
    git submodule update --init --recursive
    ./configure
    ```

### HDF5 Library Issues

If you have multiple HDF5 installations:

```bash
# Find HDF5 installations
find /usr -name "libhdf5*" 2>/dev/null
find /opt -name "libhdf5*" 2>/dev/null

# Set specific installation
export HDF5_DIR=/path/to/preferred/hdf5
./configure
```

## Next Steps

Now that cmuts is installed:

1. **[Quick Start Guide](../usage/quick-start.md)** - Learn basic usage
2. **[Command Line Options](../reference/command-line-options.md)** - Explore all available options
3. **[MPI Installation](mpi-installation.md)** - Upgrade to parallel processing (optional)
4. **[Testing](../development/testing.md)** - Run comprehensive tests

## Updating cmuts

To update to a newer version:

```bash
cd cmuts
git pull origin master
git submodule update --recursive
./configure
```

!!! note "Clean Builds"
    If you encounter issues after updating, try a clean build:
    ```bash
    rm -rf build/ bin/
    ./configure
    ```
