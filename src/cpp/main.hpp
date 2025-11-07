#ifndef _CMUTS_MAIN_HEADER_
#define _CMUTS_MAIN_HEADER_

#include "fasta.hpp"
#include "common.hpp"
#include "cmuts.hpp"
#include "hdf5.hpp"
#include "mpi.hpp"
#include "utils.hpp"
#include "options.hpp"

#ifdef MPI_BUILD
#ifdef DEBUG
const std::string PROGRAM = "cmuts core MPI (DEBUG)";
#else
const std::string PROGRAM = "cmuts core MPI";
#endif
#else
#ifdef DEBUG
const std::string PROGRAM = "cmuts core (DEBUG)";
#else
const std::string PROGRAM = "cmuts core";
#endif
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
    // Input/Output
    Arg<std::vector<std::string>> files;
    Arg<std::string> output;
    Arg<std::string> fasta;
    Arg<bool> overwrite;
    Arg<bool> rebuild;
    
    // Quality and filtering
    Arg<int> min_mapq;
    Arg<int> min_quality;
    Arg<int> min_length;
    Arg<int> max_length;
    Arg<int> max_hamming;
    Arg<float> subsample;
    
    // Processing
    Arg<int> compression;
    Arg<int> max_indel_length;
    Arg<int> chunk_size;
    Arg<int> quality_window;
    Arg<int> collapse;
    
    // Modes
    Arg<bool> joint;
    Arg<bool> tokenize;
    
    // Mutation type filters
    Arg<bool> no_mismatch;
    Arg<bool> no_insertion;
    Arg<bool> no_deletion;
    
    // Strand options
    Arg<bool> no_reverse;
    Arg<bool> only_reverse;
    
    // Deletion handling
    Arg<bool> uniform_spread;
    Arg<bool> no_spread;
    Arg<int> deletion_gap;
    Arg<bool> disable_ambiguous;
    
    // Quality filtering
    Arg<bool> no_filter_matches;
    Arg<bool> no_filter_insertions;
    Arg<bool> no_filter_deletions;

    // Base filtering
    Arg<std::string> ignore_bases;

    cmutsProgram();
};

cmutsProgram::cmutsProgram()
    : Program(PROGRAM, PROGRAM + " " + VERSION),
      files(_parser, cmuts::options::FILES.short_name, cmuts::options::FILES.long_name, cmuts::options::FILES.help, std::vector<std::string>{}),
      output(_parser, cmuts::options::OUTPUT.short_name, cmuts::options::OUTPUT.long_name, cmuts::options::OUTPUT.help),
      fasta(_parser, cmuts::options::FASTA.short_name, cmuts::options::FASTA.long_name, cmuts::options::FASTA.help),
      overwrite(_parser, cmuts::options::OVERWRITE.short_name, cmuts::options::OVERWRITE.long_name, cmuts::options::OVERWRITE.help),
      rebuild(_parser, cmuts::options::REBUILD.short_name, cmuts::options::REBUILD.long_name, cmuts::options::REBUILD.help),
      min_mapq(_parser, cmuts::options::MIN_MAPQ.short_name, cmuts::options::MIN_MAPQ.long_name, cmuts::options::MIN_MAPQ.help, cmuts::options::MIN_MAPQ.default_value, "Filtering arguments"),
      min_quality(_parser, cmuts::options::MIN_PHRED.short_name, cmuts::options::MIN_PHRED.long_name, cmuts::options::MIN_PHRED.help, cmuts::options::MIN_PHRED.default_value),
      min_length(_parser, cmuts::options::MIN_LENGTH.short_name, cmuts::options::MIN_LENGTH.long_name, cmuts::options::MIN_LENGTH.help, cmuts::options::MIN_LENGTH.default_value),
      max_length(_parser, cmuts::options::MAX_LENGTH.short_name, cmuts::options::MAX_LENGTH.long_name, cmuts::options::MAX_LENGTH.help, cmuts::options::MAX_LENGTH.default_value),
      max_hamming(_parser, cmuts::options::MAX_HAMMING.short_name, cmuts::options::MAX_HAMMING.long_name, cmuts::options::MAX_HAMMING.help, cmuts::options::MAX_HAMMING.default_value),
      subsample(_parser, cmuts::options::SUBSAMPLE.short_name, cmuts::options::SUBSAMPLE.long_name, cmuts::options::SUBSAMPLE.help, cmuts::options::SUBSAMPLE.default_value),
      compression(_parser, cmuts::options::COMPRESSION.short_name, cmuts::options::COMPRESSION.long_name, cmuts::options::COMPRESSION.help, cmuts::options::COMPRESSION.default_value),
      max_indel_length(_parser, cmuts::options::MAX_INDEL_LENGTH.short_name, cmuts::options::MAX_INDEL_LENGTH.long_name, cmuts::options::MAX_INDEL_LENGTH.help, cmuts::options::MAX_INDEL_LENGTH.default_value),
      chunk_size(_parser, cmuts::options::CHUNK_SIZE.short_name, cmuts::options::CHUNK_SIZE.long_name, cmuts::options::CHUNK_SIZE.help, cmuts::options::CHUNK_SIZE.default_value),
      quality_window(_parser, cmuts::options::QUALITY_WINDOW.short_name, cmuts::options::QUALITY_WINDOW.long_name, cmuts::options::QUALITY_WINDOW.help, cmuts::options::QUALITY_WINDOW.default_value),
      collapse(_parser, cmuts::options::COLLAPSE.short_name, cmuts::options::COLLAPSE.long_name, cmuts::options::COLLAPSE.help, cmuts::options::COLLAPSE.default_value),
      joint(_parser, cmuts::options::JOINT.short_name, cmuts::options::JOINT.long_name, cmuts::options::JOINT.help),
      tokenize(_parser, cmuts::options::TOKENIZE.short_name, cmuts::options::TOKENIZE.long_name, cmuts::options::TOKENIZE.help),
      no_mismatch(_parser, cmuts::options::NO_MISMATCH.short_name, cmuts::options::NO_MISMATCH.long_name, cmuts::options::NO_MISMATCH.help),
      no_insertion(_parser, cmuts::options::NO_INSERTION.short_name, cmuts::options::NO_INSERTION.long_name, cmuts::options::NO_INSERTION.help),
      no_deletion(_parser, cmuts::options::NO_DELETION.short_name, cmuts::options::NO_DELETION.long_name, cmuts::options::NO_DELETION.help),
      no_reverse(_parser, cmuts::options::NO_REVERSE.short_name, cmuts::options::NO_REVERSE.long_name, cmuts::options::NO_REVERSE.help),
      only_reverse(_parser, cmuts::options::ONLY_REVERSE.short_name, cmuts::options::ONLY_REVERSE.long_name, cmuts::options::ONLY_REVERSE.help),
      uniform_spread(_parser, cmuts::options::UNIFORM_SPREAD.short_name, cmuts::options::UNIFORM_SPREAD.long_name, cmuts::options::UNIFORM_SPREAD.help),
      no_spread(_parser, cmuts::options::NO_SPREAD.short_name, cmuts::options::NO_SPREAD.long_name, cmuts::options::NO_SPREAD.help),
      deletion_gap(_parser, cmuts::options::DELETION_GAP.short_name, cmuts::options::DELETION_GAP.long_name, cmuts::options::DELETION_GAP.help, cmuts::options::DELETION_GAP.default_value),
      disable_ambiguous(_parser, cmuts::options::DISABLE_AMBIGUOUS.short_name, cmuts::options::DISABLE_AMBIGUOUS.long_name, cmuts::options::DISABLE_AMBIGUOUS.help),
      no_filter_matches(_parser, cmuts::options::NO_FILTER_MATCHES.short_name, cmuts::options::NO_FILTER_MATCHES.long_name, cmuts::options::NO_FILTER_MATCHES.help),
      no_filter_insertions(_parser, cmuts::options::NO_FILTER_INSERTIONS.short_name, cmuts::options::NO_FILTER_INSERTIONS.long_name, cmuts::options::NO_FILTER_INSERTIONS.help),
      no_filter_deletions(_parser, cmuts::options::NO_FILTER_DELETIONS.short_name, cmuts::options::NO_FILTER_DELETIONS.long_name, cmuts::options::NO_FILTER_DELETIONS.help),
      ignore_bases(_parser, cmuts::options::IGNORE_BASES.short_name, cmuts::options::IGNORE_BASES.long_name, cmuts::options::IGNORE_BASES.help, cmuts::options::IGNORE_BASES.default_value)
{}

#endif
