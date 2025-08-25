#ifndef _MAIN_HEADER
#define _MAIN_HEADER

#include "fasta.hpp"
#include "common.hpp"
#include "cmuts.hpp"
#include "hdf5.hpp"
#include "mpi.hpp"
#include "utils.hpp"

#ifdef MPI_BUILD
const std::string PROGRAM = "cmuts MPI";
#else
const std::string PROGRAM = "cmuts";
#endif

static inline bool __mpi_build() {

    #ifdef MPI_BUILD
    return true;
    #else
    return false;
    #endif

}

class cmutsProgram : public Program {
public:

    Arg<std::vector<std::string>> files;
    Arg<std::string> output;
    Arg<std::string> fasta;

    Arg<bool> joint;
    Arg<bool> lowmem;
    Arg<bool> uniform_spread;
    Arg<bool> mutation_spread;

    Arg<bool> overwrite;
    Arg<int> compression;
    Arg<int> min_mapq;
    Arg<int> min_quality;
    Arg<int> max_indel_length;
    Arg<int> chunk_size;
    Arg<int> min_length;
    Arg<int> max_length;
    Arg<int> quality_window;
    Arg<int> collapse;
    Arg<bool> no_mismatch;
    Arg<bool> no_insertion;
    Arg<bool> no_deletion;
    Arg<bool> tokenize;
    Arg<float> subsample;
    Arg<bool> filter_coverage;
    Arg<bool> contiguous_ambiguous;
    Arg<bool> disable_ambiguous;

    cmutsProgram();

};

const std::string FILES_SHORT_NAME = "";
const std::string FILES_LONG_NAME = "files";
const std::string FILES_HELP = "The input SAM/BAM/CRAM files.";
const std::vector<std::string> FILES_DEFAULT;

const std::string OUTPUT_SHORT_NAME = "-o";
const std::string OUTPUT_LONG_NAME = "--output";
const std::string OUTPUT_HELP = "The output HDF5 file.";

const std::string FASTA_SHORT_NAME = "-f";
const std::string FASTA_LONG_NAME = "--fasta";
const std::string FASTA_HELP = "The reference FASTA file.";

const std::string OVERWRITE_SHORT_NAME = "";
const std::string OVERWRITE_LONG_NAME = "--overwrite";
const std::string OVERWRITE_HELP = "Overwrite an existing HDF5 file.";

const std::string COMPRESSION_SHORT_NAME = "-c";
const std::string COMPRESSION_LONG_NAME = "--compression";
const int COMPRESSION_DEFAULT = 3;
const std::string COMPRESSION_HELP = "Compression level of the HDF5 output (0-9).";

const std::string MIN_PHRED_SHORT_NAME = "";
const std::string MIN_PHRED_LONG_NAME = "--min-phred";
const int MIN_PHRED_DEFAULT = 20;
const std::string MIN_PHRED_HELP = "PHRED score threshold for base processing.";

const std::string MIN_MAPQ_SHORT_NAME = "";
const std::string MIN_MAPQ_LONG_NAME = "--min-mapq";
const int MIN_MAPQ_DEFAULT = 20;
const std::string MIN_MAPQ_HELP = "Mapping quality threshold for alignment processing.";

const std::string MAX_INDEL_LENGTH_SHORT_NAME = "";
const std::string MAX_INDEL_LENGTH_LONG_NAME = "--max-indel-length";
const int MAX_INDEL_LENGTH_DEFAULT = 10;
const std::string MAX_INDEL_LENGTH_HELP = "The longest indels to consider.";

const std::string CHUNK_SIZE_SHORT_NAME = "";
const std::string CHUNK_SIZE_LONG_NAME = "--chunk-size";
const int CHUNK_SIZE_DEFAULT = 128;
const std::string CHUNK_SIZE_HELP = "The number of references to process at a time per thread.";

const std::string MIN_LENGTH_SHORT_NAME = "";
const std::string MIN_LENGTH_LONG_NAME = "--min-length";
const int MIN_LENGTH_DEFAULT = 2;
const std::string MIN_LENGTH_HELP = "Minimum length for alignment processing.";

const std::string MAX_LENGTH_SHORT_NAME = "";
const std::string MAX_LENGTH_LONG_NAME = "--max-length";
const int MAX_LENGTH_DEFAULT = 10000;
const std::string MAX_LENGTH_HELP = "Maximum length for alignment processing.";

const std::string JOINT_SHORT_NAME = "";
const std::string JOINT_LONG_NAME = "--joint";
const std::string JOINT_HELP = "Compute the joint distribution of mutations.";

const std::string LOW_MEM_SHORT_NAME = "";
const std::string LOW_MEM_LONG_NAME = "--low-mem";
const std::string LOW_MEM_HELP = "Compute modification locations and coverage only (i.e. no modification types).";

const std::string UNIFORM_SPREAD_SHORT_NAME = "";
const std::string UNIFORM_SPREAD_LONG_NAME = "--uniform-spread";
const std::string UNIFORM_SPREAD_HELP = "Uniformly spread out ambiguous deletions.";

const std::string MUTATION_SPREAD_SHORT_NAME = "";
const std::string MUTATION_SPREAD_LONG_NAME = "--mutation-spread";
const std::string MUTATION_SPREAD_HELP = "Spread out ambiguous deletions according to the current mutation profile.";

const std::string QUALITY_WINDOW_SHORT_NAME = "";
const std::string QUALITY_WINDOW_LONG_NAME = "--quality-window";
const int QUALITY_WINDOW_DEFAULT = 2;
const std::string QUALITY_WINDOW_HELP = "Check the quality of each base in a window of this size around each base.";

const std::string COLLAPSE_SHORT_NAME = "";
const std::string COLLAPSE_LONG_NAME = "--collapse";
const int COLLAPSE_DEFAULT = 2;
const std::string COLLAPSE_HELP = "Collapse modifications within this distance of each other in a given read.";

const std::string NO_MISMATCH_SHORT_NAME = "";
const std::string NO_MISMATCH_LONG_NAME = "--no-mismatches";
const std::string NO_MISMATCH_HELP = "Do not count mismatches as modifications.";

const std::string NO_INSERTION_SHORT_NAME = "";
const std::string NO_INSERTION_LONG_NAME = "--no-insertions";
const std::string NO_INSERTION_HELP = "Do not count insertions as modifications.";

const std::string NO_DELETION_SHORT_NAME = "";
const std::string NO_DELETION_LONG_NAME = "--no-deletions";
const std::string NO_DELETION_HELP = "Do not count deletions as modifications.";

const std::string TOKENIZE_SHORT_NAME = "";
const std::string TOKENIZE_LONG_NAME = "--tokenize";
const std::string TOKENIZE_HELP = "Tokenize the reference sequences.";

const std::string SUBSAMPLE_SHORT_NAME = "";
const std::string SUBSAMPLE_LONG_NAME = "--subsample";
const float SUBSAMPLE_DEFAULT = 1.0;
const std::string SUBSAMPLE_HELP = "Randomly choose to use a read with this probability.";

const std::string FILTER_COVERAGE_SHORT_NAME = "";
const std::string FILTER_COVERAGE_LONG_NAME = "--filter-coverage";
const std::string FILTER_COVERAGE_HELP = "Apply the same filters to matches as are applied to modifications.";

const std::string CONTIGUOUS_AMBIGUOUS_SHORT_NAME = "";
const std::string CONTIGUOUS_AMBIGUOUS_LONG_NAME = "--contiguous-ambiguous";
const std::string CONTIGUOUS_AMBIGUOUS_HELP = "Allow only contiguous regions to be considered ambiguous deletions.";

const std::string DISABLE_AMBIGUOUS_SHORT_NAME = "";
const std::string DISABLE_AMBIGUOUS_LONG_NAME = "--disable-ambiguous";
const std::string DISABLE_AMBIGUOUS_HELP = "Disable the ambiguous delection detection algorithm, relying on the deletion provided by the alignment.";


cmutsProgram::cmutsProgram()
    : Program(PROGRAM, PROGRAM + " " + Utils::_get_version()),
      files(_parser, FILES_SHORT_NAME, FILES_LONG_NAME, FILES_HELP, FILES_DEFAULT),
      output(_parser, OUTPUT_SHORT_NAME, OUTPUT_LONG_NAME, OUTPUT_HELP),
      fasta(_parser, FASTA_SHORT_NAME, FASTA_LONG_NAME, FASTA_HELP),
      joint(_parser, JOINT_SHORT_NAME, JOINT_LONG_NAME, JOINT_HELP),
      lowmem(_parser, LOW_MEM_SHORT_NAME, LOW_MEM_LONG_NAME, LOW_MEM_HELP),
      uniform_spread(_parser, UNIFORM_SPREAD_SHORT_NAME, UNIFORM_SPREAD_LONG_NAME, UNIFORM_SPREAD_HELP),
      mutation_spread(_parser, MUTATION_SPREAD_SHORT_NAME, MUTATION_SPREAD_LONG_NAME, MUTATION_SPREAD_HELP),
      overwrite(_parser, OVERWRITE_SHORT_NAME, OVERWRITE_LONG_NAME, OVERWRITE_HELP),
      compression(_parser, COMPRESSION_SHORT_NAME, COMPRESSION_LONG_NAME, COMPRESSION_HELP, COMPRESSION_DEFAULT),
      min_mapq(_parser, MIN_MAPQ_SHORT_NAME, MIN_MAPQ_LONG_NAME, MIN_MAPQ_HELP, MIN_MAPQ_DEFAULT),
      min_quality(_parser, MIN_PHRED_SHORT_NAME, MIN_PHRED_LONG_NAME, MIN_PHRED_HELP, MIN_PHRED_DEFAULT),
      max_indel_length(_parser, MAX_INDEL_LENGTH_SHORT_NAME, MAX_INDEL_LENGTH_LONG_NAME, MAX_INDEL_LENGTH_HELP, MAX_INDEL_LENGTH_DEFAULT),
      chunk_size(_parser, CHUNK_SIZE_SHORT_NAME, CHUNK_SIZE_LONG_NAME, CHUNK_SIZE_HELP, CHUNK_SIZE_DEFAULT),
      min_length(_parser, MIN_LENGTH_SHORT_NAME, MIN_LENGTH_LONG_NAME, MIN_LENGTH_HELP, MIN_LENGTH_DEFAULT),
      max_length(_parser, MAX_LENGTH_SHORT_NAME, MAX_LENGTH_LONG_NAME, MAX_LENGTH_HELP, MAX_LENGTH_DEFAULT),
      quality_window(_parser, QUALITY_WINDOW_SHORT_NAME, QUALITY_WINDOW_LONG_NAME, QUALITY_WINDOW_HELP, QUALITY_WINDOW_DEFAULT),
      collapse(_parser, COLLAPSE_SHORT_NAME, COLLAPSE_LONG_NAME, COLLAPSE_HELP, COLLAPSE_DEFAULT),
      no_mismatch(_parser, NO_MISMATCH_SHORT_NAME, NO_MISMATCH_LONG_NAME, NO_MISMATCH_HELP),
      no_insertion(_parser, NO_INSERTION_SHORT_NAME, NO_INSERTION_LONG_NAME, NO_INSERTION_HELP),
      no_deletion(_parser, NO_DELETION_SHORT_NAME, NO_DELETION_LONG_NAME, NO_DELETION_HELP),
      tokenize(_parser, TOKENIZE_SHORT_NAME, TOKENIZE_LONG_NAME, TOKENIZE_HELP),
      subsample(_parser, SUBSAMPLE_SHORT_NAME, SUBSAMPLE_LONG_NAME, SUBSAMPLE_HELP, SUBSAMPLE_DEFAULT),
      filter_coverage(_parser, FILTER_COVERAGE_SHORT_NAME, FILTER_COVERAGE_LONG_NAME, FILTER_COVERAGE_HELP),
      contiguous_ambiguous(_parser, CONTIGUOUS_AMBIGUOUS_SHORT_NAME, CONTIGUOUS_AMBIGUOUS_LONG_NAME, CONTIGUOUS_AMBIGUOUS_HELP),
      disable_ambiguous(_parser, DISABLE_AMBIGUOUS_SHORT_NAME, DISABLE_AMBIGUOUS_LONG_NAME, DISABLE_AMBIGUOUS_HELP)
{}


#endif
