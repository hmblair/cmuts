```markdown
# Requirements

This page outlines the system requirements and dependencies needed to build and run cmuts.

## System Requirements

### Operating Systems

cmuts is supported on:

- **Linux** (Ubuntu 20.04+, CentOS 7+, etc.)
- **macOS** (10.15+)
- **Windows** (via WSL2 recommended)

### Hardware Requirements

=== "Minimum"
    - **CPU**: 2 cores
    - **RAM**: 4 GB
    - **Storage**: 1 GB free space

=== "Recommended"
    - **CPU**: 8+ cores (for parallel processing)
    - **RAM**: 16+ GB (for large datasets)
    - **Storage**: 10+ GB free space
    - **SSD**: For improved I/O performance

## Core Dependencies

### Required for Basic Installation

| Dependency | Minimum Version | Purpose |
|------------|----------------|---------|
| **cmake** | 3.29+ | Build system |
| **samtools** | 1.10+ | SAM/BAM/CRAM file handling |
| **HTSlib** | 1.10+ | High-throughput sequencing data |
| **HDF5** | 1.10+ | Output file format |
| **autoconf** | 2.69+ | Building htscodecs dependency |

### Additional for MPI Support

| Dependency | Purpose |
|------------|---------|
| **OpenMPI** | Parallel processing |
| **HDF5 with parallel support** | Parallel HDF5 I/O |

## Installation Methods

### Package Managers

=== "Homebrew (macOS/Linux)"
    ```bash
    brew install cmake samtools hdf5 autoconf
    
    # For MPI support
    brew install hdf5-mpi open-mpi
    ```

=== "Ubuntu/Debian"
    ```bash
    sudo apt update
    sudo apt install cmake samtools libhts-dev libhdf5-dev autoconf
    
    # For MPI support
    sudo apt install libhdf5-openmpi-dev openmpi-bin libopenmpi-dev
    ```

=== "CentOS/RHEL"
    ```bash
    sudo yum install cmake samtools htslib-devel hdf5-devel autoconf
    
    # For MPI support  
    sudo yum install hdf5-openmpi-devel openmpi-devel
    ```

### Stanford Sherlock Users

If you're running on Stanford's Sherlock cluster:

```bash
ml load hdf5/1.14.4
ml load biology samtools/1.16.1  
ml load cmake/3.31.4

# For MPI support
ml load openmpi/4.1.2
```

!!! note "Module Availability"
    Module versions may change. Use `ml avail` to see currently available versions.

## Python Dependencies (for cmuts-normalize)

The normalization component requires:

- **Python**: 3.10 or higher
- **Packages**: Listed in `requirements.txt`

```bash
pip install -r requirements.txt
```

## Environment Variables

### HDF5 Configuration

If cmake has trouble finding your HDF5 installation:

```bash
export HDF5_DIR=/path/to/hdf5/installation
```

### For Custom Installations

```bash
# Add to your shell profile (.bashrc, .zshrc, etc.)
export CMAKE_PREFIX_PATH="/usr/local:$CMAKE_PREFIX_PATH"
export PKG_CONFIG_PATH="/usr/local/lib/pkgconfig:$PKG_CONFIG_PATH"
```

## Verification

### Check Dependencies

Verify your installations:

```bash
# Check versions
cmake --version        # Should be 3.29+
samtools --version     # Should show version info
pkg-config --modversion hdf5    # Should show HDF5 version

# Check MPI (if needed)
mpirun --version       # Should show OpenMPI version
```

### Test HDF5

```bash
# Test basic HDF5 functionality
h5dump --version
```

## Common Issues

!!! warning "HDF5 Detection Problems"
    If cmake cannot find HDF5, try:
    ```bash
    export HDF5_DIR=$(brew --prefix hdf5)  # macOS with Homebrew
    # or
    export HDF5_DIR=/usr/lib/x86_64-linux-gnu/hdf5/serial  # Ubuntu
    ```

!!! tip "Samtools Version"
    Older versions of samtools may not support all CRAM features. We recommend version 1.15 or newer for best compatibility.

!!! info "MPI vs Single-threaded"
    You can install and use the single-threaded version first, then add MPI support later if needed. The basic version is sufficient for most use cases.

## Next Steps

Once you have all dependencies installed:

1. **[Basic Installation](basic-installation.md)** - Build the single-threaded version
2. **[MPI Installation](mpi-installation.md)** - Add parallel processing support
3. **[Quick Start Guide](../usage/quick-start.md)** - Begin using cmuts
