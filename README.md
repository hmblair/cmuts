## Overview

`cmuts` is a mutation-counting program that is made for wranging experimental data from modern NGS technologies. It features
* Fast, compiled C++ execution with native multithreading support
* Streamed IO and direct output to compressed HDF5 files
* Handling of arbitrary-length ambiguous deletions, including deletion spreading

## Installation

To build, you will need installed copies of
   1. `samtools` and `HTSlib`.
   2. `HDF5` with parallel support enabled.
   3. `OpenMPI`.

These must be visible to `cmake`. Stanford users: if you are running on Sherlock, then
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
   1. A FASTA file of reference sequences, as specified by the `-f` flag.
   2. One or more SAM/BAM/CRAM files of aligned reads.

Note that an index (.fai, .bai, or .crai) for each of these files will be built automatically. If the file is not sorted, `samtools sort` will be called.

In addition to the respective datasets below, the output HDF5 file will also have a `sequence` dataset containing the FASTA sequences encoded as 8-bit integers.

### Standard Mode

Counting the modifications present in a single aligned and sorted HTS file can be achieved by calling
```
cmuts -o out.h5 -f seq.fasta sorted.bam
```
To use multiple threads, invoke `mpirun -np THREADS` first:
```
mpirun -np 8 cmuts -o out.h5 -f seq.fasta sorted.bam
```
Multiple HTS files with the same reference can be processed at the same time as well:
```
mpirun -np 8 cmuts -o out.h5 -f seq.fasta sorted1.bam sorted2.bam ...
```
The output file will contain one dataset per input, with name given by the path of the file (without any extension). Each is of shape $`n \times l \times 4 \times 6`$, with dimension 2 specifying the original base and dimension 3 the mutated base. Here, deletions and insertions are considered a mutation and are in the final two rows of the array.


### Joint Mode

To compute the joint distribution of modifications over all positions, the `--joint` flag can be passed.
```
cmuts --joint -o out.h5 -f seq.fasta sorted.bam
```
The output file will contain one dataset per input with the same naming scheme as in the standard mode. It will be of shape $`n \times l \times l \times 4`$, with the final dimension specifying `$p_{0,0}, p_{0,1}, p_{1,0}$`, and `$p_{1,1}$` respectively.
