# MPI Installation

This guide covers installing cmuts with MPI (Message Passing Interface) support for parallel processing. MPI installation enables cmuts to utilize multiple CPU cores and potentially multiple nodes for faster processing of large datasets.

## When to Use MPI

Consider MPI installation if you have:

- **Large BAM/CRAM files** (>5 GB)
- **Multiple files** to process in batch
- **Multi-core systems** (8+ cores recommended)
- **High-throughput** analysis requirements
- **Cluster computing** environment available

## Prerequisites

### Basic Dependencies

All [basic installation requirements](requirements.md) plus:

| Component | Purpose |
|-----------|---------|
| **OpenMPI** | MPI implementation |
| **HDF5 with parallel support** | Parallel HDF5 I/O |

### Verify MPI Installation

```bash
# Check MPI installation
mpirun --version
mpicc --version

# Check parallel HDF5
h5pcc --version  # Should be available if HDF5 has parallel support
```

## Installation Methods

### Method 1: Package Managers (Recommended)

=== "Homebrew (macOS)"
    ```bash
    # Install parallel HDF5 and OpenMPI
    brew install hdf5-mpi open-mpi
    
    # Verify installation
    brew list hdf5-mpi
    which mpirun
    ```

=== "Ubuntu/Debian"
    ```bash
    # Install MPI and parallel HDF5
    sudo apt update
    sudo apt install openmpi-bin libopenmpi-dev \
                     libhdf5-openmpi-dev libhdf5-dev
    
    # Verify installation  
    dpkg -l | grep -E "(openmpi|hdf5)"
    ```

=== "CentOS/RHEL"
    ```bash
    # Install MPI and parallel HDF5
    sudo yum install openmpi openmpi-devel \
                     hdf5-openmpi hdf5-openmpi-devel
    
    # Load MPI module (if using modules)
    module load mpi/openmpi-x86_64
    ```

### Method 2: Stanford Sherlock

```bash
# Load required modules
ml load openmpi/4.1.2
ml load hdf5/1.14.4  # Ensure this has parallel support
ml load cmake/3.31.4
ml load biology samtools/1.16.1

# Verify parallel HDF5
ml show hdf5  # Should mention parallel or MPI support
```

## Building with MPI

### 1. Clone Repository

```bash
git clone --recurse-submodules https://github.com/hmblair/cmuts
cd cmuts
```

### 2. Configure for MPI

Run the configure script with the MPI flag:

```bash
./configure --mpi
```

This will:

- Detect MPI compiler wrappers (`mpicc`, `mpicxx`)
- Configure HDF5 with parallel support
- Build cmuts with MPI threading capabilities
- Enable parallel I/O optimizations

### 3. Verify MPI Build

```bash
# Check that MPI version was built
file bin/cmuts  # Should show it's dynamically linked
ldd bin/cmuts | grep mpi  # Should show MPI libraries

# Test MPI functionality
mpirun -np 2 bin/cmuts --help
```

## Environment Configuration

### MPI Environment Variables

Set these for optimal performance:

```bash
# Recommended MPI settings
export OMPI_MCA_btl_vader_single_copy_mechanism=none
export OMPI_MCA_mpi_warn_on_fork=0

# For large datasets
export OMPI_MCA_io_romio_cb_read=enable
export OMPI_MCA_io_romio_cb_write=enable
```

### HDF5 Parallel Configuration

```bash
# Ensure HDF5 finds MPI
export HDF5_CC=$(which mpicc)
export HDF5_CLINKER=$(which mpicc)

# For debugging HDF5 parallel I/O (optional)
export HDF5_DISABLE_VERSION_CHECK=1
```

## Testing MPI Installation

### Basic MPI Test

```bash
# Test with 4 processes
mpirun -np 4 cmuts --help

# Should output help message 4 times (one per process)
```

### Performance Test

```bash
# Run built-in performance tests
./tests/profile --mpi

# Compare single-threaded vs MPI performance
time cmuts -o test.h5 -f reference.fasta input.bam
time mpirun -np 8 cmuts -o test_mpi.h5 -f reference.fasta input.bam
```

### Comprehensive Test Suite

```bash
# Run tests with MPI
./tests/run --mpi

# Should show successful completion of all test cases
```

## Usage Examples

### Basic Parallel Execution

```bash
# Use 8 MPI processes
mpirun -np 8 cmuts -o output.h5 -f reference.fasta input.bam

# Process multiple files in parallel
mpirun -np 8 cmuts -o output.h5 -f reference.fasta \
    file1.bam file2.bam file3.bam
```

### Cluster Usage

For SLURM-based clusters:

```bash
#!/bin/bash
#SBATCH --job-name=cmuts_analysis
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=16
#SBATCH --time=04:00:00
#SBATCH --mem=64G

module load openmpi hdf5 samtools

mpirun cmuts -o results.h5 -f reference.fasta large_dataset.bam
```

## Troubleshooting

### Common MPI Issues

!!! failure "MPI not found during configure"
    **Solution**: Ensure MPI is in PATH and try explicit paths
    ```bash
    export CC=$(which mpicc)
    export CXX=$(which mpicxx)
    ./configure --mpi
    ```

!!! failure "HDF5 lacks parallel support"
    **Error**: `HDF5 was not built with parallel I/O support`
    
    **Solution**: Install HDF5 with MPI support
    ```bash
    # Homebrew
    brew uninstall hdf5
    brew install hdf5-mpi
    
    # Ubuntu
    sudo apt install libhdf5-openmpi-dev
    ```

!!! failure "Runtime MPI errors"
    **Error**: `ORTE_ERROR_LOG: Not found`
    
    **Solution**: Check MPI installation and PATH
    ```bash
    which mpirun
    mpirun --version
    echo $PATH | grep -o mpi
    ```

### Performance Issues

!!! warning "Poor scaling"
    If MPI doesn't improve performance:
    
    - **Check I/O bottlenecks**: Use SSD storage
    - **Optimize process count**: Try `nproc / 2` processes
    - **Check memory usage**: Each process needs adequate RAM
    - **Profile with**: `./tests/profile --mpi -np 4,8,16`

### HDF5 Parallel I/O Issues

!!! info "Debugging HDF5 parallel problems"
    ```bash
    # Enable HDF5 debug output
    export HDF5_DEBUG=all
    mpirun -np 4 cmuts [options]
    
    # Check HDF5 parallel capabilities
    h5pcc -showconfig | grep -i parallel
    ```

## Optimal Usage Patterns

### Process Count Guidelines

| Dataset Size | Recommended Processes |
|-------------|----------------------|
| < 1 GB | 2-4 processes |
| 1-10 GB | 4-8 processes |
| 10-50 GB | 8-16 processes |
| > 50 GB | 16-32 processes |

!!! tip "Process Count Selection"
    Start with `min(file_size_GB, num_cores)` processes and adjust based on performance.

### Memory Considerations

```bash
# Monitor memory usage
mpirun -np 8 cmuts --chunk-size 64 -o output.h5 -f ref.fasta input.bam

# For memory-constrained systems
mpirun -np 4 cmuts --low-mem --chunk-size 32 ...
```

## Performance Optimization

### Recommended Settings

```bash
# Optimized MPI execution
mpirun -np 8 \
    --map-by core \
    --bind-to core \
    cmuts -o output.h5 -f reference.fasta input.bam \
    --chunk-size 256 \
    --no-insertions \
    --filter-coverage
```

### I/O Optimization

- **Use SSD storage** for input and output files
- **Place temporary files** on local fast storage
- **Avoid network filesystems** for intensive I/O when possible

## Next Steps

With MPI cmuts installed:

1. **[Performance Testing](../development/performance.md)** - Benchmark your installation
2. **[Advanced Usage](../usage/advanced.md)** - Explore parallel processing options
3. **[Cluster Integration](../usage/cluster.md)** - Deploy on computing clusters
4. **[Optimization Guide](../usage/optimization.md)** - Tune performance parameters

## Updating MPI Installation

```bash
cd cmuts
git pull origin main
git submodule update --recursive

# Clean build for MPI
rm -rf build/ bin/
./configure --mpi
```
