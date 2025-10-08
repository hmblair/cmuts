#ifndef TESTS_HEADER
#define TESTS_HEADER

#include "cmuts.hpp"
#include "fasta.hpp"
#include "common.hpp"
#include "utils.hpp"
#include "hdf5.hpp"
#include "mpi.hpp"

class Params {
public:

    int64_t min_mapq;
    int64_t min_phred;
    int64_t min_length;
    int64_t max_length;
    int64_t max_indel_length;
    int64_t quality_window;
    int64_t collapse;
    bool mismatches;
    bool insertions;
    bool deletions;

};

const std::string PROGRAM = "cmuts generate";

class testsProgram : public Program {
public:

    Arg<int> length;
    Arg<int> queries;
    Arg<int> references;
    Arg<std::string> out_fasta;
    Arg<std::string> out_sam;
    Arg<std::string> out_h5;
    Arg<int> min_mapq;
    Arg<int> min_phred;
    Arg<int> min_length;
    Arg<int> max_length;
    Arg<int> max_indel_length;
    Arg<int> collapse;
    Arg<int> quality_window;
    Arg<bool> no_mismatch;
    Arg<bool> no_insertion;
    Arg<bool> no_deletion;

    testsProgram();

};

const std::string LENGTH_SHORT_NAME = "";
const std::string LENGTH_LONG_NAME = "--length";
const std::string LENGTH_HELP = "The length of the reference sequences.";

const std::string MIN_LENGTH_SHORT_NAME = "";
const std::string MIN_LENGTH_LONG_NAME = "--min-length";
const std::string MIN_LENGTH_HELP = "The smallest query sequences to consider when counting modifications.";
const int MIN_LENGTH_DEFAULT = 2;

const std::string MAX_LENGTH_SHORT_NAME = "";
const std::string MAX_LENGTH_LONG_NAME = "--max-length";
const std::string MAX_LENGTH_HELP = "The longest query sequences to consider when counting modifications.";

const std::string QUERIES_SHORT_NAME = "";
const std::string QUERIES_LONG_NAME = "--queries";
const std::string QUERIES_HELP = "The number of queries to generate per reference.";

const std::string REFERENCES_SHORT_NAME = "";
const std::string REFERENCES_LONG_NAME = "--references";
const std::string REFERENCES_HELP = "The number of references to generate.";

const std::string OUT_FASTA_SHORT_NAME = "";
const std::string OUT_FASTA_LONG_NAME = "--out-fasta";
const std::string OUT_FASTA_HELP = "The file to store the references in.";

const std::string OUT_SAM_SHORT_NAME = "";
const std::string OUT_SAM_LONG_NAME = "--out-sam";
const std::string OUT_SAM_HELP = "The file to store the queries in.";

const std::string OUT_H5_SHORT_NAME = "";
const std::string OUT_H5_LONG_NAME = "--out-h5";
const std::string OUT_H5_HELP = "The file to store the expected modifications in.";

const std::string MIN_MAPQ_SHORT_NAME = "";
const std::string MIN_MAPQ_LONG_NAME = "--min-mapq";
const std::string MIN_MAPQ_HELP = "The minimum quality to consider a read.";

const std::string MIN_PHRED_SHORT_NAME = "";
const std::string MIN_PHRED_LONG_NAME = "--min-phred";
const std::string MIN_PHRED_HELP = "The minimum quality to consider a base.";


const std::string MAX_INDEL_SHORT_NAME = "";
const std::string MAX_INDEL_LONG_NAME = "--max-indel-length";
const std::string MAX_INDEL_HELP = "Skip indels longer than this.";

const std::string WINDOW_SHORT_NAME = "";
const std::string WINDOW_LONG_NAME = "--quality-window";
const std::string WINDOW_HELP = "The number of neighbouring bases to consider when calculating PHRED scores.";

const std::string COLLAPSE_SHORT_NAME = "";
const std::string COLLAPSE_LONG_NAME = "--collapse";
const std::string COLLAPSE_HELP = "The minimum number of bases between modifications to consider them consecutive.";
const int COLLAPSE_DEFAULT = 2;

const std::string NO_DEL_SHORT_NAME = "";
const std::string NO_DEL_LONG_NAME = "--no-deletions";
const std::string NO_DEL_HELP = "Do not count deletions as modifications.";

const std::string NO_INS_SHORT_NAME = "";
const std::string NO_INS_LONG_NAME = "--no-insertions";
const std::string NO_INS_HELP = "Do not count insertions as modifications.";

const std::string NO_MIS_SHORT_NAME = "";
const std::string NO_MIS_LONG_NAME = "--no-mismatches";
const std::string NO_MIS_HELP = "Do not count mismatches as modifications.";




testsProgram::testsProgram()
    : Program(PROGRAM, VERSION),
      length(_parser, LENGTH_SHORT_NAME, LENGTH_LONG_NAME, LENGTH_HELP),
      queries(_parser, QUERIES_SHORT_NAME, QUERIES_LONG_NAME, QUERIES_HELP),
      references(_parser, REFERENCES_SHORT_NAME, REFERENCES_LONG_NAME, REFERENCES_HELP),
      out_fasta(_parser, OUT_FASTA_SHORT_NAME, OUT_FASTA_LONG_NAME, OUT_FASTA_HELP),
      out_sam(_parser, OUT_SAM_SHORT_NAME, OUT_SAM_LONG_NAME, OUT_SAM_HELP),
      out_h5(_parser, OUT_H5_SHORT_NAME, OUT_H5_LONG_NAME, OUT_H5_HELP),
      min_mapq(_parser, MIN_MAPQ_SHORT_NAME, MIN_MAPQ_LONG_NAME, MIN_MAPQ_HELP),
      min_phred(_parser, MIN_PHRED_SHORT_NAME, MIN_PHRED_LONG_NAME, MIN_PHRED_HELP),
      min_length(_parser, MIN_LENGTH_SHORT_NAME, MIN_LENGTH_LONG_NAME, MIN_LENGTH_HELP, MIN_LENGTH_DEFAULT),
      max_length(_parser, MAX_LENGTH_SHORT_NAME, MAX_LENGTH_LONG_NAME, MAX_LENGTH_HELP),
      max_indel_length(_parser, MAX_INDEL_SHORT_NAME, MAX_INDEL_LONG_NAME, MAX_INDEL_HELP),
      collapse(_parser, COLLAPSE_SHORT_NAME, COLLAPSE_LONG_NAME, COLLAPSE_HELP, COLLAPSE_DEFAULT),
      quality_window(_parser, WINDOW_SHORT_NAME, WINDOW_LONG_NAME, WINDOW_HELP),
      no_mismatch(_parser, NO_MIS_SHORT_NAME, NO_MIS_LONG_NAME, NO_MIS_HELP),
      no_insertion(_parser, NO_INS_SHORT_NAME, NO_INS_LONG_NAME, NO_INS_HELP),
      no_deletion(_parser, NO_DEL_SHORT_NAME, NO_DEL_LONG_NAME, NO_DEL_HELP) {}

#endif
