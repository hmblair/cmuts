#ifndef _CMUTS_MAIN_HEADER_
#define _CMUTS_MAIN_HEADER_

#include "fasta.hpp"
#include "common.hpp"
#include "cmuts.hpp"
#include "hdf5.hpp"
#include "mpi.hpp"
#include "utils.hpp"
#include "generated/cmuts_args.hpp"

static inline bool __mpi_build() {
    #ifdef MPI_BUILD
    return true;
    #else
    return false;
    #endif
}

class cmutsProgram : public Program {
public:
    CMUTSPROGRAM_ARG_MEMBERS

    cmutsProgram()
        : Program(CMUTSPROGRAM_PROGRAM_NAME, CMUTSPROGRAM_PROGRAM_VERSION),
          CMUTSPROGRAM_ARG_INIT
    {}
};

#endif
