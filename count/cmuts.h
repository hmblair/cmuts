#include <stdint.h>
#include <string.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
// getopt
#include <getopt.h>
#include <errno.h>
// For printing integers with commas
#include <locale.h>
// Processing indexed BAM and FASTA files
#include "htsutils.h"
// HDF5
#include "h5utils.h"
// Binary sequences
#include "bseq.h"

#define VERSION "0.2.0"
// The name of the datasets in the HDF5 file
#define MUTATION_DS "mutations"
#define INSERTION_DS "insertions"
// The number of dimensions of the respective
// datasets
#define NUM_MUTDIMS 4
#define NUM_INSDIMS 3
// The default out file
#define DEFAULT_OUT "count.h5"
// For printing in colour
#define RED "\033[91m"
#define END_CLR "\033[0m"

// For MPI
#define ROOT_PROC 0

