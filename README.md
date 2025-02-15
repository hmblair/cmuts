## Overview

`cmuts` is a mutation-counting program that is made for wranging experimental data from modern NGS technologies. It features
* Fast, compiled C execution with native multithreading support
* Streamed IO and direct output to compressed HDF5 files
* Handling of arbitrary-length ambiguous indels, including indel spreading

## Installation

To build, you will need installed copies of
   1. `samtools` and `HTSlib`.
   2. `HDF5` with parallel support enabled.
   3. `OpenMPI`.

These must be visible to `cmake`. If you are running on Sherlock, then 
```
ml load hdf5/1.14.4
ml load biology samtools/1.16.1
ml load cmake/3.24.2
```
is sufficient.

Installation is done by cloning the repository and running the installation script.
```
git clone https://github.com/hmblair/cmuts
cd cmuts
./configure
```
Don't forget to add the `./bin` directory to your path.

If `cmake` has issues finding the HDF5 installation or you want to use a specific one, set the `HDF5_DIR` environment variable to the desired installation directory and it will be used instead.

## Usage

To run, you will need:
   1. An indexed FASTA file, as specified by the `-f` flag. This means that there is a corresponding .fai file, which can be generated with `samtools faidx`.
   2. An indexed SAM/BAM/CRAM file. This means that there is a corresponding .bai file, which can be generated by first sorting the file and then calling `samtools index`.

See `cmuts --help` for more runtime options.

The output file will contain two datasets, called "mutations" and "insertions". The former is of shape $`n \times l \times 4 \times 5`$, with dimension 2 specifying the original base and dimension 3 the mutated base. Here, deletions are considered a mutation and are in the final row of the array. The latter is of shape $`n \times l \times 4`$, with the final dimension denoting the type of the inserted base.
