#ifndef TESTS_HEADER
#define TESTS_HEADER

#include "hts.hpp"
#include "cmuts.hpp"
#include "utils.hpp"
#include "hdf5.hpp"
#include "mpi.hpp"

const std::string PROGRAM = "cmuts-test";
const std::string VERSION = "1.0.0";

class testsProgram : public Program {
public:

    Arg<int> length;
    Arg<int> queries;
    Arg<int> references;
    Arg<std::string> out_fasta;
    Arg<std::string> out_sam;
    Arg<std::string> out_h5;
    Arg<int> min_quality;
    Arg<int> min_length;
    Arg<int> window;

    testsProgram();

};

// Variables for "--length" argument
const std::string LENGTH_SHORT_NAME = "";
const std::string LENGTH_LONG_NAME = "--length";
const std::string LENGTH_HELP = "The length of the reference sequences.";

const std::string MIN_LENGTH_SHORT_NAME = "";
const std::string MIN_LENGTH_LONG_NAME = "--min-length";
const std::string MIN_LENGTH_HELP = "The smallest query sequences to consider when counting modifications.";
const int MIN_LENGTH_DEFAULT = 2;

// Variables for "--queries" argument
const std::string QUERIES_SHORT_NAME = "";
const std::string QUERIES_LONG_NAME = "--queries";
const std::string QUERIES_HELP = "The number of queries to generate per reference.";

// Variables for "--references" argument
const std::string REFERENCES_SHORT_NAME = "";
const std::string REFERENCES_LONG_NAME = "--references";
const std::string REFERENCES_HELP = "The number of references to generate.";

// Variables for "--out-fasta" argument
const std::string OUT_FASTA_SHORT_NAME = "";
const std::string OUT_FASTA_LONG_NAME = "--out-fasta";
const std::string OUT_FASTA_HELP = "The file to store the references in.";

// Variables for "--out-sam" argument
const std::string OUT_SAM_SHORT_NAME = "";
const std::string OUT_SAM_LONG_NAME = "--out-sam";
const std::string OUT_SAM_HELP = "The file to store the queries in.";

// Variables for "--out-h5" argument
const std::string OUT_H5_SHORT_NAME = "";
const std::string OUT_H5_LONG_NAME = "--out-h5";
const std::string OUT_H5_HELP = "The file to store the expected modifications in.";

// Variables for "--min-quality" argument
const std::string MIN_QUALITY_SHORT_NAME = "";
const std::string MIN_QUALITY_LONG_NAME = "--min-quality";
const std::string MIN_QUALITY_HELP = "The minimum quality to consider a base or read.";

// Variables for "--quality-window" argument
const std::string WINDOW_SHORT_NAME = "";
const std::string WINDOW_LONG_NAME = "--quality-window";
const std::string WINDOW_HELP = "The number of neighbouring bases to consider when calculating PHRED scores.";


testsProgram::testsProgram()
    : Program(PROGRAM, VERSION),
      length(_parser, LENGTH_SHORT_NAME, LENGTH_LONG_NAME, LENGTH_HELP),
      queries(_parser, QUERIES_SHORT_NAME, QUERIES_LONG_NAME, QUERIES_HELP),
      references(_parser, REFERENCES_SHORT_NAME, REFERENCES_LONG_NAME, REFERENCES_HELP),
      out_fasta(_parser, OUT_FASTA_SHORT_NAME, OUT_FASTA_LONG_NAME, OUT_FASTA_HELP),
      out_sam(_parser, OUT_SAM_SHORT_NAME, OUT_SAM_LONG_NAME, OUT_SAM_HELP),
      out_h5(_parser, OUT_H5_SHORT_NAME, OUT_H5_LONG_NAME, OUT_H5_HELP),
      min_quality(_parser, MIN_QUALITY_SHORT_NAME, MIN_QUALITY_LONG_NAME, MIN_QUALITY_HELP),
      min_length(_parser, MIN_LENGTH_SHORT_NAME, MIN_LENGTH_LONG_NAME, MIN_LENGTH_HELP, MIN_LENGTH_DEFAULT),
      window(_parser, WINDOW_SHORT_NAME, WINDOW_LONG_NAME, WINDOW_HELP) {}

#endif
