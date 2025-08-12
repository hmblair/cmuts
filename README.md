## Overview

`cmuts` is a program for counting mutations and computing reactivity profiles in MaP-seq experiments. It features
* Fast, compiled C++ code with native multithreading support
* Streamed IO and direct output to compressed HDF5 files
* Handling of arbitrary-length ambiguous deletions, including mutation-informed deletion spreading


# Installation

## cmuts

To build single-threaded `cmuts`, you will need installed copies of
   1. `cmake >= 3.29`.
   2. `samtools` and `HTSlib`
   3. `HDF5`

Building `htscodecs`, used for reading CRAM files, also requires `autoreconf`.

The easiest way to install these dependencies is via `brew`, where they are available under `cmake`, `samtools`, `hdf5`, and `autoconf` respectively.

If `cmake` has issues finding the HDF5 installation or you want to use a specific one, set the `HDF5_DIR` environment variable to the desired installation directory and it will be used instead.

Stanford users: if you are running on Sherlock, then
```
ml load hdf5/1.14.4
ml load biology samtools/1.16.1
ml load cmake/3.31.4
```
is sufficient.

Installation is done by cloning the repository and running the installation script.
```
git clone --recurse-submodules https://github.com/hmblair/cmuts
cd cmuts
./configure
```
This will download and build the `htscodecs` dependency as well. Don't forget to add the `./bin` directory to your path.


## cmuts MPI

Running `./configure --mpi` will build the parallel version of `cmuts`. This requires `OpenMPI` and a copy of `HDF5` with parallel support enabled. This is available under `brew` as `hdf5-mpi`.


# Usage

To run, you will need:
   1. A FASTA file of reference sequences, as specified by the `-f` option.
   2. One or more SAM/BAM/CRAM files of aligned reads.

Note that an index (.bai or .crai) for the SAM/BAM/CRAM file will be created. If the file is not sorted, `samtools sort` will be called to facilitate this. A custom binary .binfa file will be created from the FASTA file.


## Modification Counting

Counting the modifications present in a single aligned HTS file can be achieved by calling
```
cmuts -o out.h5 -f seq.fasta sorted.bam
```
If built with MPI, invoke `mpirun -np THREADS` first to use multiple threads:
```
mpirun -np 8 cmuts -o out.h5 -f seq.fasta sorted.bam
```
Multiple SAM/BAM/CRAM files with the same reference can be processed at the same time as well:
```
mpirun -np 8 cmuts -o out.h5 -f seq.fasta sorted1.bam sorted2.bam ...
```
The output file will contain one dataset per input, with name given by the path of the file (without any extension). Each is of shape $`n \times l \times 4 \times 6`$ for $`n`$ sequences with a maximum length of $`l`$. Dimension 2 specifies the original base and dimension 3 the mutated base; here, deletions are considered a mutation and are in the second to final row of the array. The final row contains all insertions, with dimension 2 specifying the base inserted.

The following lists all additional commands available:

`--overwrite`: Overwrite any existing HDF5 file.

`--compression`: Compression level of the HDF5 output (0-9). Default: 3

`--chunk-size`: The size of the internal buffer (in references) per thread. Default: 128

`--tokenize`: Tokenize the reference sequences. See below.

`--joint`: Compute the joint distribution of mutations. See below.

`--low-mem`: Compute modification locations and coverage only (i.e. modification type is not recorded). Decreases memory usage by at least 2x.

`--min-mapq`: Mapping quality threshold for alignment processing. Default: 20

`--min-phred`: PHRED score threshold for base processing. Default: 20

`--max-indel-length`: The longest indels to consider as modifications. Default: 10

`--min-length`: Minimum length for alignment processing. Default: 2

`--max-length`: Maximum length for alignment processing. Default: 10000

`--quality-window`: Check the quality of each base in a window of this size around each base. Default: 2

`--collapse`: Collapse modifications within this distance of each other in a given read. Default: 2

`--uniform-spread`: Uniformly spread out ambiguous deletions.

`--mutation-spread`: Spread out ambiguous deletions according to the existing mutation profile.

`--no-mismatches`: Do not count mismatches as modifications.

`--no-insertions`: Do not count insertions as modifications.

`--no-deletions`: Do not count deletions as modifications.

`--subsample`: Randomly accept each read with this probability, to simulate lower-read experiments. Default: 1.0

`--filter-coverage`: Apply the same filters to matches as are applied to modifications.


## Normalization

The program `cmuts-normalize` will produce normalized reactivity profiles from the output of `cmuts`. It relies on the Python dependencies in `requirements.txt` and can be run as
```
cmuts-normalize -o reactivity.h5 --mod-ds MODS [--nomod-ds NOMODS] --out-groups OUTPUT_GROUPS INPUT.h5
```
The `--mod-ds` and `--nomod-ds` flags specify which dataset(s) in the input HDF5 file to process. The latter is optional, but if given, the number of datasets for both must match. `--out-groups` specifies which group(s) in the output HDF5 file to place the two datasets which are computed -- `reactivity` and `reads`.

Additional flags which may be useful are:

`--overwrite`: Overwrite an existing HDF5 file.

`--clip-reactivity`: Clip the reactivity values to the range $`[0,1]`$.

`--5p-primer-length`: The length of the 5' primer, which will be zeroed out. Default: 26

`--3p-primer-length`: The length of the 3' primer, which will be zeroed out. Default: 20


## Tokenization

If the `--tokenize` flag is provided, `cmuts` will tokenize the reference sequences and place them in the `sequence` dataset of the output HDF5 file. This can be done even in the absence of a SAM/BAM/CRAM file, i.e.
```
cmuts --tokenize -f seq.fasta -o out.h5
```
will fill the `sequence` dataset in `out.h5`.


## Joint Modification Counting

To compute the joint distribution of modifications over all pairs of positions, the `--joint` flag can be passed.
```
cmuts --joint -o out.h5 -f seq.fasta sorted.bam
```
The output file will contain one dataset per input with the same naming scheme as in the standard mode. It will be of shape $`n \times l \times l \times 2 \times 2`$, with the final two dimensions specifying whether each pair of positions has a modification ($`i=1`$) or no modification ($`i=0`$).


## Tests

`cmuts` has tests for the basic mutation counting features, which are build automatically with the main program. They can be run with `./tests/run`, which selects a random set of parameters to generate test cases with and then runs `cmuts`. The position and type of each match, mismatch, insertion and deletion are scored against their expected positions. Pass the `--cram` flag to test with CRAM rather than BAM files.

There is also a command to assess the performance of `cmuts`, `./tests/profile`. This requires `rf-count` and, on Mac OS, `gtime` to be installed.
