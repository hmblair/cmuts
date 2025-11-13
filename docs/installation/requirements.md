
## Dependencies

Building and running the `cmuts` pipeline requires the following packages:

 - `python>=3.9` with all packages in `requirements.txt`
 - `cmake>=3.29` with `pkg-config`
 - `autoconf`
 - `samtools` and `htslib`
 - `hdf5`
 - `omp`
 - `fastp`

The MPI build also requires

 - `openmpi`
 - `hdf5-mpi`

In addition, demultiplexing with `cmuts align` requires

 - `ultraplex`

to be installed. Note that this requires `python==3.9`, and will not install for more recent versions of `python`.

The following should install all (save for `python` and its packages) on a personal device. For installation on a managed cluster, consult the respective guidelines.

=== "Mac OS"
    ```bash
    # Install brew from https://brew.sh
    brew install cmake autoconf samtools
    # For non-MPI builds only:
    brew install hdf5
    # For MPI builds only:
    brew install openmpi hdf5-mpi
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
