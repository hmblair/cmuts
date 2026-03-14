
## Dependencies

Building and running the `cmuts` pipeline requires the following packages:

**Required:**

 - `python>=3.9` with all packages in `requirements.txt`
 - `cmake>=3.29` with `pkg-config`
 - `autoconf`, `automake`, `libtool`
 - `samtools` and `htslib`
 - `hdf5`
 - `fastp`

**Optional:**

 - `omp` — enables multithreaded pairwise counting (builds fine without it)

The MPI build also requires:

 - `openmpi`
 - `hdf5-mpi`

Demultiplexing with `cmuts align` requires:

 - `ultraplex`

Note that `ultraplex` requires `python==3.9`, and will not install for more recent versions of `python`.

The following should install all (save for `python` and its packages) on a personal device. For installation on a managed cluster, consult the respective guidelines.

=== "Mac OS"
    ```bash
    # Install brew from https://brew.sh
    # Required:
    brew install cmake autoconf automake libtool samtools
    brew install hdf5       # for non-MPI builds
    # Optional:
    brew install libomp     # enables multithreaded pairwise counting
    # For MPI builds only:
    brew install openmpi hdf5-mpi
    ```

=== "Linux (Ubuntu/Debian)"
    ```bash
    sudo apt-get update
    # Required:
    sudo apt-get install -y cmake autoconf pkg-config
    sudo apt-get install -y samtools libhts-dev
    sudo apt-get install -y libhdf5-dev
    # Optional:
    sudo apt-get install -y libomp-dev   # multithreaded pairwise counting
    # For MPI builds only:
    sudo apt-get install -y libopenmpi-dev openmpi-bin
    sudo apt-get install -y libhdf5-openmpi-dev
    ```

=== "Linux (Fedora/RHEL)"
    ```bash
    # Required:
    sudo dnf install -y cmake autoconf pkg-config
    sudo dnf install -y samtools htslib-devel
    sudo dnf install -y hdf5-devel
    # Optional:
    sudo dnf install -y libomp-devel     # multithreaded pairwise counting
    # For MPI builds only:
    sudo dnf install -y openmpi-devel hdf5-openmpi-devel
    # Load MPI module (may be required):
    module load mpi/openmpi-x86_64
    ```

!!! tip "Verify Dependencies"
    Run these commands to verify the dependencies are successfully installed:
    ```bash
    cmake --version
    pkg-config --modversion hdf5
    autoreconf --version
    samtools --version
    h5ls --version
    ```
    For MPI builds, also run:
    ```bash
    mpirun --version
    ```

### Python Dependencies

To create a `conda` environment with all required dependencies, you may run

```bash
conda create -n cmuts python=3.9
conda activate cmuts
pip3 install -r requirements.txt
conda install -c bioconda fastp
```

For the `cmuts align` dependencies required for demultiplexing, also run

```bash
conda install -c bioconda ultraplex
```

### HDF5 Configuration

If `cmake` has trouble finding your HDF5 installation, you can set

```bash
export HDF5_DIR=/path/to/hdf5/installation
```

If installed via brew, the command `brew info hdf5` may be helpful for finding the desired path.

On Linux, HDF5 is typically installed to `/usr/lib/x86_64-linux-gnu/hdf5` (Debian/Ubuntu) or `/usr/lib64` (Fedora/RHEL). You can find it with:

```bash
pkg-config --variable=libdir hdf5
```
