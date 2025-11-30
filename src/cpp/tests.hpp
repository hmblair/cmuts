#ifndef TESTS_HEADER
#define TESTS_HEADER

#include "cmuts.hpp"
#include "fasta.hpp"
#include "common.hpp"
#include "utils.hpp"
#include "hdf5.hpp"
#include "mpi.hpp"
#include "generated/tests_args.hpp"

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

class testsProgram : public Program {
public:
    TESTSPROGRAM_ARG_MEMBERS

    testsProgram()
        : Program(TESTSPROGRAM_PROGRAM_NAME, TESTSPROGRAM_PROGRAM_VERSION),
          TESTSPROGRAM_ARG_INIT
    {}
};

#endif
