#include "cmuts.h"
#include "h5utils.h"
#include "htsutils.h"
#include <complex.h>
#include <stdbool.h>

char *removeExtension(const char *filename) {

    char *newFilename = NULL;
    newFilename = strdup(filename);

    char *dot = strrchr(newFilename, '.');  // Find the last '.'
    if (dot != NULL) {
        *dot = '\0';  // Truncate the string at the last '.'
    }

    return newFilename;

}

char *replaceExtension(
    const char *filename,
    const char *ext
) {

    char *newFilename = NULL;
    newFilename = strdup(filename);

    char *dot = strrchr(newFilename, '.');  // Find the last '.'
    if (dot != NULL) {
        *dot = '\0';  // Truncate the string at the last '.'
    }
    strcat(newFilename, ext);

    return newFilename;

}

void replaceChar(char *str, char x, char y) {
    for (int i = 0; str[i] != '\0'; i++) {
        if (str[i] == x) {
            str[i] = y;
        }
    }
}

void fileCheck(
    const char *filename,
    bool overwrite,
    bool *outFileExists
) {
    // The file only exists if we find it and do
    // not overwrite it
    *outFileExists = false;
    if (!filename) {
        return;
    }
    if (access(filename, F_OK) == 0) {
        if (overwrite) {
            int status = remove(filename);
            if (status != 0) {
                perror("Error overwriting file:");
                exit(EXIT_FAILURE);
            }
        } else {
            *outFileExists = true;
        }
    }

}

void fileCheckParallel(
    const char *filename,
    int rank,
    bool overwrite,
    bool *outFileExists
) {

    // Use the root process to check if the
    // file exists
    if (rank == ROOT_PROC) {
        fileCheck(
            filename,
            overwrite,
            outFileExists
        );
    }
    // Broadcast the value to all other processes
    MPI_Bcast(
        outFileExists,
        1, MPI_C_BOOL,
        ROOT_PROC,
        MPI_COMM_WORLD
    );
    MPI_Barrier(MPI_COMM_WORLD);

}

char pluralise(int n) {
    if (n == 1) {
        return '\0';
    } else {
        return 's';
    }
}

typedef struct {
    double startTime;
    double currentTime;
} Timer;

Timer startTimer() {
    Timer timer;
    timer.startTime = MPI_Wtime();
    timer.currentTime = timer.startTime;
    return timer;
}

double endTimer(Timer timer) {
    timer.currentTime = MPI_Wtime();
    return timer.currentTime - timer.startTime;
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
    uint64_t maxReferenceLength;
} Progress;

Progress initProgress() {
    Progress progress;
    progress.numReferences = 0;
    progress.numAlignments = 0;
    progress.skippedAlignments = 0;
    progress.unmappedReads = 0;
    progress.totalAlignments = 0;
    progress.totalReferences = 0;
    progress.maxReferenceLength = 0;
    return progress;
}

static inline void printProgress(Progress progress, Timer timer) {

    float progress_f = (float)progress.numAlignments / progress.totalAlignments * 100;
    float skipped_f = (float)progress.skippedAlignments / progress.totalAlignments * 100;
    double time = endTimer(timer);
    char *timestring = getTimeString(time);

    printf("\033[A\033[A\033[A\033[F\033[A\r        Alignments processed: %.1f%%.\n", progress_f);
    printf("        Alignments skipped:   %.1f%%\n", skipped_f);
    printf("\x1b[2K");
    printf("        Time elapsed:         %s\n", timestring);
    printf("      ───────────────────────────────────────\n\n");
    // Flush to print now, and free the string
    fflush(stdout);
    free(timestring);

}

typedef struct {
    uint64_t ix;
    uint64_t size;
} BatchInfo;

BatchInfo getBatchInfo(uint64_t ix, uint64_t size, uint64_t max) {
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
    for (uint64_t ix = batch.ix; ix < batch.ix + batch.size; ix++) {

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


MPI_Comm splitProcs(
    int rank,
    int maxrank
) {

    MPI_Comm new_comm;
    if (rank < maxrank) {
        MPI_Comm_split(MPI_COMM_WORLD, 0, rank, &new_comm); // Group 0 for the first N processes
    } else {
        MPI_Comm_split(MPI_COMM_WORLD, MPI_UNDEFINED, rank, &new_comm); // Exclude others
    }

    return new_comm;

}


void printVersion() {
    printf("\n  cmuts version %s\n", VERSION);
}
const char *USAGE_STR = "Usage: mpirun -np PROCS cmuts -f FASTA [options] [BAM FILES]";
void printHelp(int STATUS) {
    printVersion();
    printf("  %s\n\n", USAGE_STR);
    printf("    --fasta=S, -f:           The input FASTA file, possibly in .gz format, and indexed by samtools faidx.\n");
    printf("    --output=S, -o:          The output HDF5 file.\n");
    printf("    --overwite:              Whether to delete the existing file, should it exist.\n");
    printf("    --compression=N, -c:     The compression level of the HDF5 file, from 0 to 9. Defaults to 3.\n");
    printf("    --spread, -s:  Spread ambiguous deletions across the region of ambiguity.\n");
    printf("    --collapse=N:            Collapse mutations which are within N bases towards the 3' end.\n");
    printf("    --min-mapping-quality=N: Skip alignments with a mapping quality below N.\n");
    printf("    --min-mut-base-quality=N:    Ignore mutations with a PHRED score below N.\n");
    printf("    --min-cov-base-quality=N:    Ignore coverage with a PHRED score below N.\n");
    printf("    --min-query-length=N:    Skip alignments with a query length below N.\n");
    printf("    --max-indel-length=N:    Ignore indels with a length longer than N.\n");
    printf("    --num-neighbours=N:      Check the PHRED score of N neighbours on each side of the base.\n");
    printf("\n");
    exit(STATUS);
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


int main(int argc, char **argv) {

    if (argc == 1) {
        printHelp(EXIT_FAILURE);
    }

    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const char *inFile = NULL;
    const char *fastaFile = NULL;
    const char *outFile = DEFAULT_OUT;

    // Flags to be passed to the dataspace constructor
    DataspaceFlags dFlags = defaultDataspaceFlags();
    // Flags to be passed to the accumulateMutations function
    CountingFlags cFlags = defaultCountingFlags();

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
        {"spread", no_argument, 0, 's'},
        {"collapse", required_argument, 0, 0},
        {"min-query-length", required_argument, 0, 0},
        {"min-mapping-quality", required_argument, 0, 0},
        {"min-mut-base-quality", required_argument, 0, 0},
        {"min-cov-base-quality", required_argument, 0, 0},
        {"max-indel-length", required_argument, 0, 0},
        {"num-neighbours", required_argument, 0, 0},
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
            if (strcmp(longOptions[optionIndex].name, "min-mut-base-quality") == 0) {
                cFlags.minBaseMutQ = safeStringToInt(optarg, "min-mut-base-quality");
            }
            if (strcmp(longOptions[optionIndex].name, "min-cov-base-quality") == 0) {
                cFlags.minBaseCovQ = safeStringToInt(optarg, "min-cov-base-quality");
            }
            if (strcmp(longOptions[optionIndex].name, "max-indel-length") == 0) {
                cFlags.maxDelLength = safeStringToInt(optarg, "max-indel-length");
            }
            if (strcmp(longOptions[optionIndex].name, "min-mapping-quality") == 0) {
                cFlags.minSeqMapQ = safeStringToInt(optarg, "min-mapping-quality");
            }
            if (strcmp(longOptions[optionIndex].name, "min-query-length") == 0) {
                cFlags.minQueryLength= safeStringToInt(optarg, "min-query-length");
            }
            if (strcmp(longOptions[optionIndex].name, "num-neighbours") == 0) {
                cFlags.numNeighboursToCheck = safeStringToInt(optarg, "num-neighbours");
            }
            break;
          case 'h':
            printHelp(EXIT_SUCCESS);
          case 'f':
            fastaFile = optarg;
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

    bool fastaFileExists;
    fileCheck(
        fastaFile,
        false,
        &fastaFileExists
    );
    if (!fastaFileExists) {
        fprintf(stderr, "No FASTA file was provided. please provide one with the -f flag.\n");
    }

    bool cFlagsValid = verifyCountingFlags(cFlags);
    if (!cFlagsValid) {
        exit(EXIT_FAILURE);
    }

    // Find out how many files we are going to be processing
    int nfiles = argc - optind;
    if (nfiles <= 0) {
        printf("No SAM/BAM/CRAM file was provided.\n");
        exit(EXIT_FAILURE);
    }

    // Check if the output file exists; if it does, check if the
    // user wants it deleted.
    bool outFileExists;
    fileCheckParallel(
        outFile,
        rank,
        overwrite,
        &outFileExists
    );

    // Create a Progress object to keep track of the
    // current progress
    Progress progress = initProgress();
    // Loop over all specified files
    for (int nfile = 0; nfile < nfiles; nfile++) {

        inFile = argv[optind + nfile];
        bool inFileExists;
        fileCheck(
            inFile,
            false,
            &inFileExists
        );
        if (!inFileExists) {
            fprintf(stderr, "The provided SAM/BAM/CRAM file \"%s\" cannot be opened.\n", inFile);
            exit(EXIT_FAILURE);
        }

        // Open an IndexedBAM handle to get some metadata
        // on the file
        IndexedBAM meta = openIndexedBAM(
            inFile,
            fastaFile,
            false
        );

        if (nfile == 0) {
            progress.totalReferences += meta.numReferences;
            progress.maxReferenceLength = meta.maxReferenceLength;
        }
        progress.totalAlignments += meta.numAlignments;
        progress.unmappedReads += meta.numUnmappedReads;

        closeIndexedBAM(meta);

    }


    setlocale(LC_NUMERIC, "");
    if (rank == ROOT_PROC) {

    printf("\n        cmuts version %s\n", VERSION);
    printf("      ───────────────────────────────────────\n");
    printf("        Total references:     %'" PRIu64 "\n", progress.totalReferences);
    printf("        Total alignments:     %'" PRIu64 "\n", progress.totalAlignments);
    printf("        Total unmapped reads: %'" PRIu64 "\n", progress.unmappedReads);
    printf("      ───────────────────────────────────────\n\n\n\n\n\n");

    }

    // Start a timer to track the time taken
    Timer timer = startTimer();

    // Either make or open the filespace using the
    // getFilespace function
    Filespace filespace = getFilespace(
        outFile,
        outFileExists,
        H5F_ACC_RDWR,
        MPI_COMM_WORLD,
        MPI_INFO_NULL
    );

    uint64_t maxIX = (progress.totalReferences / (size * CHUNK_SIZE) ) * size * CHUNK_SIZE;
    uint64_t overflowProcs = (maxIX + CHUNK_SIZE - 1) / CHUNK_SIZE;
    // MPI_Comm MPI_COMM_SPLIT = splitProcs(
    //    rank,
    //    overflowProcs
    //);
    // Get a second filespace to handle writing
    // any overflow that not all threads can
    // participate in
    // Filespace overflowFilespace = getFilespace(
    //     outFile,
    //    true,
    //    H5F_ACC_RDWR,
    //    MPI_COMM_SPLIT,
    //    MPI_INFO_NULL
    //);

    for (int nfile = 0; nfile < nfiles; nfile++) {

    inFile = argv[optind + nfile];

    // Remove the extension from the file, and replace all backslashes with
    // hyphens, to get the group name
    dFlags.group = removeExtension(inFile);
    replaceChar(dFlags.group, '/', '-');

    // Specify the size of the datasets in the HDF5 file
    hsize_t mutationDimensions[4] = {
        progress.totalReferences, 
        progress.maxReferenceLength, 
        N_BASES,
        N_DELBASES
    };
    hsize_t insertionDimensions[3] = {
        progress.totalReferences, 
        progress.maxReferenceLength, 
        N_BASES
    };
    uint64_t totalRefs = progress.totalReferences;

    if (groupExists(filespace, dFlags.group)) {
        printf("The group \"%s\" in the file \"%s\" already exists. Please use the --overwrite flag to overwrite the file, or choose another group name to write to the existing file with the --group option.\n", dFlags.group, outFile);
        closeFilespace(filespace);
        exit(EXIT_FAILURE);
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

    // Each thread gets their own local memspace to
    // write to
    Memspace mutationMemspace = initMemspace(mutationDataspace, CHUNK_SIZE);
    Memspace insertionMemspace = initMemspace(insertionDataspace, CHUNK_SIZE);

    // Each thread gets their own file handle,
    // so that data may be loaded in parallel
    IndexedBAM ixBAM = openIndexedBAM(
        inFile,
        fastaFile,
        true
    );
    IndexedFASTA ixFASTA = openIndexedFASTA(fastaFile);
    // Set the reference genome for the BAM file

    // The main loop
    for (
        uint64_t refIX = rank * CHUNK_SIZE;
        refIX < maxIX;
        refIX += size * CHUNK_SIZE
    ) {

        BatchInfo batch = getBatchInfo(
            refIX,
            CHUNK_SIZE,
            maxIX
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

        uint64_t localSum = 0;
        MPI_Reduce(
            &count.numAlignments, 
            &localSum, 
            1, MPI_INT, MPI_SUM, 
            0, MPI_COMM_WORLD
        );
        uint64_t localSkipped = 0;
        MPI_Reduce(
            &count.numSkipped, 
            &localSkipped, 
            1, MPI_INT, MPI_SUM, 
            0, MPI_COMM_WORLD
        );

        if (rank == ROOT_PROC) {
            progress.numAlignments += localSum;
            progress.skippedAlignments += localSkipped;
            printProgress(progress, timer);
        }

    }

    // Close this thread's BAM and FASTA handles
    closeIndexedBAM(ixBAM);
    closeIndexedFASTA(ixFASTA);
    // Free this thread's memspaces
    freeMemspace(mutationMemspace);
    freeMemspace(insertionMemspace);
    // Close the remaining HDF5 objects
    closeDataspace(mutationDataspace);
    closeDataspace(insertionDataspace);


    // TODO: process the final loop of references
    // if not divisible by CHUNK_SIZE

    }

    closeFilespace(filespace);
    // closeFilespace(overflowFilespace);
    H5close();

    MPI_Finalize();

    return EXIT_SUCCESS;

}
