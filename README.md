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

# Usage

To run, you will need:
   1. A FASTA file of reference sequences, as specified by the `-f` flag.
   2. One or more SAM/BAM/CRAM files of aligned reads.

Note that an index (.fai, .bai, or .crai) for each of these files will be built automatically. If the file is not sorted, `samtools sort` will be called.

In addition to the respective datasets below, the output HDF5 file will also have a `sequence` dataset containing the FASTA sequences encoded as 8-bit integers.

## Modification Counting

Counting the modifications present in a single aligned HTS file can be achieved by calling
```
cmuts -o out.h5 -f seq.fasta sorted.bam
```
To use multiple threads, invoke `mpirun -np THREADS` first:
```
mpirun -np 8 cmuts -o out.h5 -f seq.fasta sorted.bam
```
Multiple SAM/BAM/CRAM files with the same reference can be processed at the same time as well:
```
mpirun -np 8 cmuts -o out.h5 -f seq.fasta sorted1.bam sorted2.bam ...
```
The output file will contain one dataset per input, with name given by the path of the file (without any extension). Each is of shape $`n \times l \times 4 \times 6`$ for $`n`$ sequences with a maximum length of $`l`$. Dimension 2 specifies the original base and dimension 3 the mutated base; here, deletions are considered a mutation and are in the second to final row of the array. The final row contains all insertions, with dimension 2 specifying the base inserted.

The following lists all additional commands available:

`--overwrite`: Overwrite an existing HDF5 file.

`--compression`: Compression level of the HDF5 output (0-9). Default: 3

`--min-phred`: PHRED score threshold for base processing. Default: 20

`--min-mapq`: Mapping quality threshold for alignment processing. Default: 20

`--max-indel-length`: The longest indels to consider. Default: 10

`--chunk-size`: The number of references to process at a time per thread. Default: 128

`--min-length`: Minimum length for alignment processing. Default: 2

`--max-length`: Maximum length for alignment processing. Default: 10000

`--quality-window`: Check the quality of each base in a window of this size around each base. Default: 2

`--joint`: Compute the joint distribution of mutations. See below.

`--fast`: Compute modification locations and coverage only.

`--very-fast`: Compute modification locations only.

`--spread`: Spread out ambiguous deletions.

## Joint Mutation Counting

To compute the joint distribution of modifications over all positions, the `--joint` flag can be passed.
```
cmuts --joint -o out.h5 -f seq.fasta sorted.bam
```
The output file will contain one dataset per input with the same naming scheme as in the standard mode. It will be of shape $`n \times l \times l \times 2 \times 2`$, with the final two dimensions specifying whether each position has a modification ($`i=1`$) or no modification ($`i=0`$).

## Normalization

The program `cmuts-normalize` will produce normalized reactivity profiles from the output of `cmuts`. It relies on the dependencies in `requiresments.txt` and can be run as
```
cmuts-normalize -o reactivity.h5 --mod-ds MODS [--nomod-ds NOMODS] --out-group OUTPUT_GROUPS INPUT.h5
```
The `--mod-ds` and `--nomod-ds` flags specify which datasets in the input HDF5 file to process. The latter is optional, but if given, the number of datasets for both must match. `--out-group` specifies which group in the output HDF5 file to place the two datasets, `reactivity` and `reads`; the same rules apply for it too.

Additional flags which may be useful are:

`overwrite`: Overwrite an existing HDF5 file.

`--clip-reactivity`: Clip the reactivity values to the range $`[0,1]`$.

`--5p-primer-length`: The length of the 5' primer, which will be zeroed out. Default: 26

`--3p-primer-length`: The length of the 3' primer, which will be zeroed out. Default: 20

## Tests

`cmuts` has tests for the basic mutation counting features, which are build automatically with the main program. They can be run with `./tests/accuracy/run`, which selects a random set of parameters to generate test cases with and then runs `cmuts`. The position and type of each match, mismatch, insertion and deletion are scored.
