#include <stdint.h>
#include <string.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
// OMP
#include <omp.h>
// getopt
#include <getopt.h>
#include <errno.h>
// For printing integers with commas
#include <locale.h>
// Processing indexed BAM and FASTA files
#include "htsutils.h"
// HDF5
#include "h5utils.h"
// Constants
#include "cmuts.h"

#define MIN(a,b) ((a) < (b) ? (a) : (b))

bool endsWith(const char *str, const char *suffix) {

    size_t stlen = strlen(str);
    size_t sulen = strlen(suffix);
    if (stlen < sulen) {
        return false;
    }
    for (int i = 0; i < sulen; i++) {
        if (str[i + stlen - sulen] != suffix[i]) {
            return false;
        }
    }
    return true;

}

char pluralise(int n) {
    if (n == 1) {
        return '\0';
    } else {
        return 's';
    }
}

typedef struct {
    struct timespec startTime;
    struct timespec endTime;
} Timer;

Timer startTimer() {
    Timer timer;
    clock_gettime(CLOCK_MONOTONIC, &timer.startTime);
    return timer;
}

double endTimer(Timer timer) {
    clock_gettime(CLOCK_MONOTONIC, &timer.endTime);
    return (timer.endTime.tv_sec - timer.startTime.tv_sec) + (timer.endTime.tv_nsec - timer.startTime.tv_nsec) / 1e9;
}

#define MINUTE 60
#define HOUR 3600
char *getTimeString(double time) {
    char *str;
    if (time < MINUTE) {
        asprintf(&str, "%.2f seconds.", time);
        return str;
    } else if (time < HOUR) {
        asprintf(&str, "%.2f minutes.", time / MINUTE);
        return str;
    } else {
        asprintf(&str, "%.2f hours.", time / HOUR);
        return str;
    }
}

typedef struct {
    uint64_t numReferences;
    uint64_t numAlignments;
    uint64_t skippedAlignments;
    uint64_t totalReferences;
    uint64_t totalAlignments;
    uint64_t unmappedReads;
} Progress;

Progress initProgress() {
    Progress progress;
    progress.numReferences = 0;
    progress.numAlignments = 0;
    progress.skippedAlignments = 0;
    return progress;
}

static inline void printProgress(Progress progress, Timer timer) {

    float progress_f = (float)progress.numAlignments / progress.totalAlignments * 100;
    float skipped_f = (float)progress.skippedAlignments / progress.totalAlignments * 100;
    double time = endTimer(timer);
    char *timestring = getTimeString(time);

    #pragma omp critical 
    {
        printf("\033[A\033[A\033[A\033[F\033[A\r        Alignments processed: %.1f%%.\n", progress_f);
        printf("        Alignments skipped:   %.1f%%\n", skipped_f);
        printf("\x1b[2K");
        printf("        Time elapsed:         %s\n", timestring);
        printf("      ───────────────────────────────────────\n\n");
        fflush(stdout);
    }

}

typedef struct {
    int ix;
    int size;
} BatchInfo;

BatchInfo getBatchInfo(int ix, int size, int max) {
    BatchInfo batch;
    batch.ix = ix;
    batch.size = MIN(size, max - ix);
    return batch;
}

AlignmentCount batchedCountMutations(
    IndexedBAM ixBAM,
    IndexedFASTA ixFASTA,
    BatchInfo batch,
    Memspace mutationMemspace,
    Memspace insertionMemspace,
    CountingFlags cFlags
    ) {

    // Resize the memspaces to align with the size of the 
    // amount of data to be processed
    resizeMemspace(&mutationMemspace, batch.size);
    resizeMemspace(&insertionMemspace, batch.size);

    // Zero out the slabs
    fillMemspace(mutationMemspace, 0.0f);
    fillMemspace(insertionMemspace, 0.0f);

    // Get a copy of the pointers to the beginning of the mutation
    // and insertion arrays.
    float *mutations = mutationMemspace.data;
    float *insertions = insertionMemspace.data;

    AlignmentCount count;
    count.numAlignments = 0;
    count.numSkipped = 0;
    for (int32_t ix = batch.ix; ix < batch.ix + batch.size; ix++) {

        // Get the reference sequence for this index
        char *sequence = getSequence(ixFASTA, ix);
        // Accumulate mutations in the mutation and insertion
        // arrays, and get the number of alignments processed
        AlignmentCount count_tmp = countMutations(
            mutations,
            insertions,
            sequence,
            ixBAM,
            ix,
            cFlags
        );
        count.numAlignments += count_tmp.numAlignments;
        count.numSkipped += count_tmp.numSkipped;
        // Free the reference sequence
        free(sequence);

        // Get the length of the reference sequence
        int referenceLength = ixBAM.header->target_len[ix];
        // Incrememt the pointers to point to the next reference
        // sequence
        mutations += referenceLength * N_BASES * N_DELBASES;
        insertions += referenceLength * N_BASES;

    }

    return count;

}

void printVersion() {
    printf("\n  cmuts version %s\n", VERSION);
}
const char *USAGE_STR = "Usage: cmuts [options] in.bam";
void printHelp(int STATUS) {
    printVersion();
    printf("  %s\n\n", USAGE_STR);
    printf("    --fasta=S, -f:           The input FASTA file, possibly in .gz format, and indexed by samtools faidx.\n");
    printf("    --output=S, -o:          The output HDF5 file.\n");
    printf("    --group=S:               Put the output datasets inside this group in the HDF5 file. This allows for appending to existing HDF5 files.\n");
    printf("    --overwite:              Whether to delete the existing file, should it exist.\n");
    printf("    --compression=N, -c:     The compression level of the HDF5 file, from 0 to 9.\n");
    printf("    --num-threads=N, -t:     The number of threads to use.\n");
    printf("    --spread, -s:  Spread ambiguous deletions across the region of ambiguity.\n");
    printf("    --collapse=N:            Collapse mutations which are within N bases towards the 3' end.\n");
    printf("    --min-mapping-quality=N: Skip alignments with a mapping quality below N.\n");
    printf("    --min-base-quality=N:    Ignore mutations with a PHRED score below N.\n");
    printf("    --min-query-length=N:    Skip alignments with a query length below N.\n");
    printf("    --max-indel-length=N:    Ignore indels with a length longer than N.\n");
    printf("    --num-neighbours=N:      Check the PHRED score of N neighbours on each side of the base.\n");
    printf("\n");
    exit(STATUS);
}

bool checkAccess(const char *filename) {
    if (access(filename, F_OK) == 0) {
        return true;
    }
    return false;
}

int safeStringToInt(char *str, const char *msg) {

    char *endptr;
    errno = 0;
    int num = strtol(str, &endptr, 10);
    if (errno == ERANGE || *endptr != '\0' || str == endptr) {
        printf("\"%s\" is not a valid input to %s.\n", str, msg);
        exit(EXIT_FAILURE);
    }

    return num;

}

void verifyFile(const char *filename, const char *type) {
    if (!filename) {
        fprintf(stderr, "A %s file must be specified with the -f flag.\n", type);
        exit(EXIT_FAILURE);
    } else if (!checkAccess(filename)) {
        fprintf(stderr, "There was an error opening the %s file \"%s\".\n", type, filename);
        exit(EXIT_FAILURE);
    }
}

void verifyThreads(int n) {
    if (n < 1) {
        printf("At least one thread must be used; %d were specified.\n", n);
        exit(EXIT_FAILURE);
    }
}

int main(int argc, char **argv) {

    if (argc == 1) {
        printHelp(EXIT_FAILURE);
    }

    const char *inFile = NULL;
    const char *fastaFile = NULL;
    const char *outFile = DEFAULT_OUT;

    // Flags to be passed to the dataspace constructor
    DataspaceFlags dFlags = defaultDataspaceFlags();
    // Flags to be passed to the accumulateMutations function
    CountingFlags cFlags = defaultCountingFlags();

    int NUM_THREADS = 1;
    bool overwrite = false;

    char c;
    extern char *optarg; extern int optind;
    int optionIndex = 0;
    static const struct option longOptions[] = {
        {"help", no_argument, NULL, 'h'},
        {"version", no_argument, NULL, 0},
        {"fasta", required_argument, NULL, 'f'},
        {"output", required_argument, NULL, 'o'},
        {"compression", required_argument, NULL, 'c'},
        {"overwrite", no_argument, 0, 0},
        {"num-threads", required_argument, 0, 't'},
        {"spread", no_argument, 0, 's'},
        {"collapse", required_argument, 0, 0},
        {"min-query-length", required_argument, 0, 0},
        {"min-mapping-quality", required_argument, 0, 0},
        {"min-base-quality", required_argument, 0, 0},
        {"max-indel-length", required_argument, 0, 0},
        {"num-neighbours", required_argument, 0, 0},
        {"group", required_argument, 0, 0},
        {0, 0, 0, 0}
    };

    while ((c = getopt_long(argc, argv, "hsf:o:t:c:", longOptions, &optionIndex)) != -1) {
        switch (c) {
          case 0:
            if (strcmp(longOptions[optionIndex].name, "version") == 0) {
                printVersion();
                printf("\n");
                exit(EXIT_SUCCESS);
            }
            if (strcmp(longOptions[optionIndex].name, "collapse") == 0) {
                cFlags.collapseMutations = safeStringToInt(optarg, "collapse");
            }
            if (strcmp(longOptions[optionIndex].name, "overwrite") == 0) {
                overwrite = true;
            }
            if (strcmp(longOptions[optionIndex].name, "min-base-quality") == 0) {
                cFlags.minBaseQ = safeStringToInt(optarg, "min-base-quality");
            }
            if (strcmp(longOptions[optionIndex].name, "max-indel-length") == 0) {
                cFlags.maxDelLength = safeStringToInt(optarg, "max-indel-length");
            }
            if (strcmp(longOptions[optionIndex].name, "min-mapping-quality") == 0) {
                cFlags.minMapQ = safeStringToInt(optarg, "min-mapping-quality");
            }
            if (strcmp(longOptions[optionIndex].name, "min-query-length") == 0) {
                cFlags.minQueryLength= safeStringToInt(optarg, "min-query-length");
            }
            if (strcmp(longOptions[optionIndex].name, "num-neighbours") == 0) {
                cFlags.numNeighboursToCheck = safeStringToInt(optarg, "num-neighbours");
            }
            if (strcmp(longOptions[optionIndex].name, "group") == 0) {
                dFlags.group = optarg;
            }
            break;
          case 'h':
            printHelp(EXIT_SUCCESS);
          case 'f':
            fastaFile = optarg;
            break;
          case 't':
            NUM_THREADS = safeStringToInt(optarg, "num-threads");
            break;
          case 'o':
            outFile = optarg;
            break;
          case 's':
            cFlags.spreadDeletions = true;
            break;
          case 'c':
            dFlags.compression = safeStringToInt(optarg, "compression");
            break;
          default:
            exit(EXIT_FAILURE);
          }
    }

    if (optind == argc) {
        printf("No SAM/BAM/CRAM file was provided.\n");
        exit(EXIT_FAILURE);
    } else if (optind < argc - 1) {
        printHelp(EXIT_FAILURE);
    } else {
        inFile = argv[optind];
        if (!checkAccess(inFile)) {
            fprintf(stderr, "The provided SAM/BAM/CRAM file \"%s\" cannot be opened.\n", inFile);
            exit(EXIT_FAILURE);
        }
    }

    verifyFile(fastaFile, "FASTA");
    verifyThreads(NUM_THREADS);
    bool cFlagsValid = verifyCountingFlags(cFlags);
    if (!cFlagsValid) {
        exit(EXIT_FAILURE);
    }

    // Check that the outfile is an HDF5 file
    if (!endsWith(outFile, HDF5_SUFFIX)) {
        fprintf(stderr, "The output file must be an HDF5 (%s) file.\n", HDF5_SUFFIX);
        exit(EXIT_FAILURE);
    }
    // Check if the output file exists; if it does, check if the
    // user wants it deleted
    bool outFileExists = false;
    if (checkAccess(outFile)) {
        if (overwrite) {
            int status = remove(outFile);
            if (status != 0) {
                printf("There was an error deleting the existing output file \"%s\".\n", outFile);
                exit(EXIT_FAILURE);
            }
        } else {
            // If the group is not specified, we cannot proceed
            if (!dFlags.group) {
                printf("The file \"%s\" already exists, and a group has not been specified. Please use the --overwrite flag to overwrite the file, or specify a group to write to the existing file with the --group option.\n", outFile);
                exit(EXIT_FAILURE);
            } 
            // Else, we check if the specified group exists in the file
            else {
                // Open the filespace
                Filespace filespace = openFilespace(outFile, H5F_ACC_RDONLY);
                if (groupExists(filespace, dFlags.group)) {
                    printf("The group \"%s\" in the file \"%s\" already exists. Please use the --overwrite flag to overwrite the file, or choose another group name to write to the existing file with the --group option.\n", dFlags.group, outFile);
                    closeFilespace(filespace);
                    exit(EXIT_FAILURE);
                } 
                // Else, we need not do anything
                closeFilespace(filespace);
                outFileExists = true;
            }
        }
    }

    // Open an IndexedBAM handle to get some metadata
    // on the file
    IndexedBAM meta = openIndexedBAM(inFile);
    // Specify the size of the datasets in the HDF5 file
    hsize_t mutationDimensions[4] = {
        meta.numReferences, 
        meta.maxReferenceLength, 
        N_BASES,
        N_DELBASES
    };
    hsize_t insertionDimensions[3] = {
        meta.numReferences, 
        meta.maxReferenceLength, 
        N_BASES
    };
    // Create a Progress object to keep track of the
    // current progress
    Progress progress = initProgress();
    progress.totalReferences = meta.numReferences;
    progress.totalAlignments = meta.numAlignments;
    progress.unmappedReads = meta.numUnmappedReads;
    // Close the metadata handle
    closeIndexedBAM(meta);

    // Make an HDF5 file to store the output, or open it
    // if the file exists already.
    Filespace filespace;
    if (outFileExists) {
        filespace = openFilespace(outFile, H5F_ACC_RDWR);
    } else {
        filespace = makeFilespace(outFile);
    }

    // Make a dataset within the file to store the mutations
    Dataspace mutationDataspace = makeDataspace(
        filespace, 
        MUTATION_DS, 
        NUM_MUTDIMS, 
        mutationDimensions, 
        dFlags
    );
    // Make a dataset within the file to store the insertions
    Dataspace insertionDataspace = makeDataspace(
        filespace, 
        INSERTION_DS, 
        NUM_INSDIMS, 
        insertionDimensions, 
        dFlags
    );

    setlocale(LC_NUMERIC, "");
    printf("\n        cmuts version %s\n", VERSION);
    printf("      ───────────────────────────────────────\n");
    printf("        Total references:     %'" PRIu64 "\n", progress.totalReferences);
    printf("        Total alignments:     %'" PRIu64 "\n", progress.totalAlignments);
    printf("        Total unmapped reads: %'" PRIu64 "\n", progress.unmappedReads);
    printf("      ───────────────────────────────────────\n\n\n\n\n\n");

    // Start a timer to track the time taken
    Timer timer = startTimer();
    #pragma omp parallel num_threads(NUM_THREADS)
    {

        // Each thread gets their own local memspace to
        // write to
        Memspace mutationMemspace = initMemspace(mutationDataspace, CHUNK_SIZE);
        Memspace insertionMemspace = initMemspace(insertionDataspace, CHUNK_SIZE);

        // Each thread gets their own file handle,
        // so that data may be loaded in parallel
        IndexedBAM ixBAM = openIndexedBAM(inFile);
        IndexedFASTA ixFASTA = openIndexedFASTA(fastaFile);

        // The main loop
        #pragma omp for schedule(dynamic)
        for (int32_t refIX = 0; refIX < progress.totalReferences; refIX+=CHUNK_SIZE) {

            BatchInfo batch = getBatchInfo(
                refIX,
                CHUNK_SIZE,
                progress.totalReferences
            );

            AlignmentCount count = batchedCountMutations(
                ixBAM, 
                ixFASTA, 
                batch, 
                mutationMemspace,
                insertionMemspace,
                cFlags
            );

            writeSlab(
                mutationDataspace, 
                mutationMemspace, 
                refIX
            );
            writeSlab(
                insertionDataspace, 
                insertionMemspace, 
                refIX
            );

            #pragma omp atomic
            progress.numReferences += batch.size;
            #pragma omp atomic
            progress.numAlignments += count.numAlignments;
            #pragma omp atomic
            progress.skippedAlignments += count.numSkipped;

            printProgress(progress, timer);

        }

        // Close this thread's BAM and FASTA handles
        closeIndexedBAM(ixBAM);
        closeIndexedFASTA(ixFASTA);
        // Free this thread's memspaces
        freeMemspace(mutationMemspace);
        freeMemspace(insertionMemspace);

    }

    // Close the remaining HDF5 objects
    closeDataspace(mutationDataspace);
    closeDataspace(insertionDataspace);
    closeFilespace(filespace);
    H5close();

    return EXIT_SUCCESS;

}
