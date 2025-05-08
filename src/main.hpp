#ifndef MAIN_HEADER
#define MAIN_HEADER

#include "cmuts.hpp"
#include "hdf5.hpp"
#include "hts.hpp"
#include "mpi.hpp"
#include "utils.hpp"
#include <system_error>

const std::string PROGRAM = "cmuts";
const std::string VERSION = "1.0.0";

class cmutsProgram : public Program {
public:

    Arg<std::vector<std::string>> files;
    Arg<std::string> output;
    Arg<std::string> fasta;
    Arg<bool> overwrite;
    Arg<int> compression;
    Arg<int> min_mapq;
    Arg<int> min_quality;
    Arg<int> max_indel_length;
    Arg<int> chunk_size;
    Arg<int> min_length;
    Arg<int> max_length;
    Arg<bool> joint;
    Arg<bool> fast;
    Arg<bool> spread;
    Arg<int> quality_window;
    Arg<bool> no_mismatch;
    Arg<bool> no_insertion;
    Arg<bool> no_deletion;

    cmutsProgram();

};

// Variables for "files" argument
const std::string FILES_SHORT_NAME = "";
const std::string FILES_LONG_NAME = "files";
const std::string FILES_HELP = "The input SAM/BAM/CRAM files.";
const std::vector<std::string> FILES_DEFAULT;

// Variables for "-o", "--output" argument
const std::string OUTPUT_SHORT_NAME = "-o";
const std::string OUTPUT_LONG_NAME = "--output";
const std::string OUTPUT_HELP = "The output HDF5 file.";

// Variables for "-f", "--fasta" argument
const std::string FASTA_SHORT_NAME = "-f";
const std::string FASTA_LONG_NAME = "--fasta";
const std::string FASTA_HELP = "The reference FASTA file.";

// Variables for "--overwrite" argument
const std::string OVERWRITE_SHORT_NAME = "";
const std::string OVERWRITE_LONG_NAME = "--overwrite";
const std::string OVERWRITE_HELP = "Overwrite an existing HDF5 file.";

// Variables for "-c", "--compression" argument
const std::string COMPRESSION_SHORT_NAME = "-c";
const std::string COMPRESSION_LONG_NAME = "--compression";
const int COMPRESSION_DEFAULT = 3;
const std::string COMPRESSION_HELP = "Compression level of the HDF5 output (0-9).";

// Variables for "--min-phred" argument
const std::string MIN_PHRED_SHORT_NAME = "";
const std::string MIN_PHRED_LONG_NAME = "--min-phred";
const int MIN_PHRED_DEFAULT = 20;
const std::string MIN_PHRED_HELP = "PHRED score threshold for base processing.";

const std::string MIN_MAPQ_SHORT_NAME = "";
const std::string MIN_MAPQ_LONG_NAME = "--min-mapq";
const int MIN_MAPQ_DEFAULT = 20;
const std::string MIN_MAPQ_HELP = "Mapping quality threshold for alignment processing.";

// Variables for "--max-indel-length" argument
const std::string MAX_INDEL_LENGTH_SHORT_NAME = "";
const std::string MAX_INDEL_LENGTH_LONG_NAME = "--max-indel-length";
const int MAX_INDEL_LENGTH_DEFAULT = 10;
const std::string MAX_INDEL_LENGTH_HELP = "The longest indels to consider.";

// Variables for "--chunk-size" argument
const std::string CHUNK_SIZE_SHORT_NAME = "";
const std::string CHUNK_SIZE_LONG_NAME = "--chunk-size";
const int CHUNK_SIZE_DEFAULT = 128;
const std::string CHUNK_SIZE_HELP = "The number of references to process at a time per thread.";

// Variables for "--min-length" argument
const std::string MIN_LENGTH_SHORT_NAME = "";
const std::string MIN_LENGTH_LONG_NAME = "--min-length";
const int MIN_LENGTH_DEFAULT = 2;
const std::string MIN_LENGTH_HELP = "Minimum length for alignment processing.";

// Variables for "--max-length" argument
const std::string MAX_LENGTH_SHORT_NAME = "";
const std::string MAX_LENGTH_LONG_NAME = "--max-length";
const int MAX_LENGTH_DEFAULT = 10000;
const std::string MAX_LENGTH_HELP = "Maximum length for alignment processing.";

// Variables for "--joint" argument
const std::string JOINT_SHORT_NAME = "";
const std::string JOINT_LONG_NAME = "--joint";
const std::string JOINT_HELP = "Compute the joint distribution of mutations.";

// Variables for "--fast" argument
const std::string FAST_SHORT_NAME = "";
const std::string FAST_LONG_NAME = "--fast";
const std::string FAST_HELP = "Compute modification locations and coverage only.";

// Variables for "--fast" argument
const std::string VERY_FAST_SHORT_NAME = "";
const std::string VERY_FAST_LONG_NAME = "--very-fast";
const std::string VERY_FAST_HELP = "Compute modification locations only.";

// Variables for "--spread" argument
const std::string SPREAD_SHORT_NAME = "";
const std::string SPREAD_LONG_NAME = "--spread";
const std::string SPREAD_HELP = "Spread out ambiguous deletions.";

// Variables for "--quality-window" argument
const std::string QUALITY_WINDOW_SHORT_NAME = "";
const std::string QUALITY_WINDOW_LONG_NAME = "--quality-window";
const int QUALITY_WINDOW_DEFAULT = 2;
const std::string QUALITY_WINDOW_HELP = "Check the quality of each base in a window of this size around each base.";

const std::string NO_MISMATCH_SHORT_NAME = "";
const std::string NO_MISMATCH_LONG_NAME = "--no-mismatches";
const std::string NO_MISMATCH_HELP = "Do not count mismatches as modifications.";

const std::string NO_DELETION_SHORT_NAME = "";
const std::string NO_DELETION_LONG_NAME = "--no-insertions";
const std::string NO_DELETION_HELP = "Do not count mismatches as modifications.";

const std::string NO_INSERTION_SHORT_NAME = "";
const std::string NO_INSERTION_LONG_NAME = "--no-deletions";
const std::string NO_INSERTION_HELP = "Do not count mismatches as modifications.";

cmutsProgram::cmutsProgram()
    : Program(PROGRAM, VERSION),
      files(_parser, FILES_SHORT_NAME, FILES_LONG_NAME, FILES_HELP, FILES_DEFAULT),
      output(_parser, OUTPUT_SHORT_NAME, OUTPUT_LONG_NAME, OUTPUT_HELP),
      fasta(_parser, FASTA_SHORT_NAME, FASTA_LONG_NAME, FASTA_HELP),
      overwrite(_parser, OVERWRITE_SHORT_NAME, OVERWRITE_LONG_NAME, OVERWRITE_HELP),
      compression(_parser, COMPRESSION_SHORT_NAME, COMPRESSION_LONG_NAME, COMPRESSION_HELP, COMPRESSION_DEFAULT),
      min_mapq(_parser, MIN_MAPQ_SHORT_NAME, MIN_MAPQ_LONG_NAME, MIN_MAPQ_HELP, MIN_MAPQ_DEFAULT),
      min_quality(_parser, MIN_PHRED_SHORT_NAME, MIN_PHRED_LONG_NAME, MIN_PHRED_HELP, MIN_PHRED_DEFAULT),
      max_indel_length(_parser, MAX_INDEL_LENGTH_SHORT_NAME, MAX_INDEL_LENGTH_LONG_NAME, MAX_INDEL_LENGTH_HELP, MAX_INDEL_LENGTH_DEFAULT),
      chunk_size(_parser, CHUNK_SIZE_SHORT_NAME, CHUNK_SIZE_LONG_NAME, CHUNK_SIZE_HELP, CHUNK_SIZE_DEFAULT),
      min_length(_parser, MIN_LENGTH_SHORT_NAME, MIN_LENGTH_LONG_NAME, MIN_LENGTH_HELP, MIN_LENGTH_DEFAULT),
      max_length(_parser, MAX_LENGTH_SHORT_NAME, MAX_LENGTH_LONG_NAME, MAX_LENGTH_HELP, MAX_LENGTH_DEFAULT),
      joint(_parser, JOINT_SHORT_NAME, JOINT_LONG_NAME, JOINT_HELP),
      fast(_parser, FAST_SHORT_NAME, FAST_LONG_NAME, FAST_HELP),
      spread(_parser, SPREAD_SHORT_NAME, SPREAD_LONG_NAME, SPREAD_HELP),
      quality_window(_parser, QUALITY_WINDOW_SHORT_NAME, QUALITY_WINDOW_LONG_NAME, QUALITY_WINDOW_HELP, QUALITY_WINDOW_DEFAULT),
      no_mismatch(_parser, NO_MISMATCH_SHORT_NAME, NO_MISMATCH_LONG_NAME, NO_MISMATCH_HELP),
      no_deletion(_parser, NO_DELETION_SHORT_NAME, NO_DELETION_LONG_NAME, NO_DELETION_HELP),
      no_insertion(_parser, NO_INSERTION_SHORT_NAME, NO_INSERTION_LONG_NAME, NO_INSERTION_HELP)
{}


#endif
