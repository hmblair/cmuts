#ifndef CMUTS_HEADER
#define CMUTS_HEADER

#include "tiny.hpp"
#include "hdf5.hpp"
#include "mpi.hpp"
#include "utils.hpp"
#include <random>

const int64_t N_BASES    = 4;
const int64_t N_DELBASES = 6;

const std::string MODIFICATION_DS = "modifications";
const std::string JOINT_DS        = "joint";

const hts_pos_t NOMOD = 0;
const hts_pos_t MOD   = 1;

const hts_pos_t MOD_MOD     = 0;
const hts_pos_t MOD_NOMOD   = 1;
const hts_pos_t NOMOD_MOD   = 2;
const hts_pos_t NOMOD_NOMOD = 3;

namespace cmuts {





//
// Enums
//





enum class Mode {

    // Count mutation types, locations, and coverage
    // Size: (N, L, 4, 6)
    Normal,
    // Count mutation locations and coverage only
    // Size: (N, L, 2)
    LowMem,
    // Compute the joint distribution of modifications
    // Size: (N, L, L, 4)
    Joint,
    // Tokenize the reference sequences
    // Size: (N, L)
    Tokenize

};


constexpr size_t _ndims(Mode mode) {

    switch (mode) {
        case Mode::Normal:   { return 4; }
        case Mode::LowMem:   { return 3; }
        case Mode::Joint:    { return 5; }
        case Mode::Tokenize: { return 2; }
    }

}


enum class Spread{

    None,
    Uniform,
    MutationInformed

};






//
// Stats
//





class Stats {
private:

    int64_t _processed  = 0;
    int64_t _skipped    = 0;
    int64_t _aligned    = 0;
    int64_t _unaligned  = 0;
    int64_t _references = 0;

    const MPI::Manager& _mpi;

    Utils::Line _print_processed = Utils::Line("Reads processed", "%");
    Utils::Line _print_skipped   = Utils::Line("Reads skipped", "%");
    Utils::Line _print_elapsed   = Utils::Line("Time elapsed");

public:

    Stats(
        int64_t aligned,
        int64_t unaligned,
        int64_t references,
        const MPI::Manager& mpi
    );

    void processed();
    void skipped();
    void aggregate();
    void header() const;
    void body() const;

};





//
// Params
//





class Params {
public:

    Mode mode;
    Spread spread;
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
    float subsample;

};





//
// Main
//





class __Main {
protected:

    TinyHTS::File& file;
    TinyHTS::FASTA& fasta;
    HDF5::File& hdf5;
    const MPI::Manager& mpi;
    const Params& params;
    MPI::Chunk chunk;
    Stats& stats;

public:

    __Main(
        TinyHTS::File& file,
        TinyHTS::FASTA& fasta,
        HDF5::File& hdf5,
        const MPI::Manager& mpi,
        const Params& params,
        Stats& stats
    );

    virtual ~__Main() = default;

    virtual void run() = 0;

};

template <typename dtype, Mode mode, Spread spread>
class Main : public __Main {
private:

    HDF5::Memspace<dtype, _ndims(mode)> memspace;

public:

    Main(
        TinyHTS::File& file,
        TinyHTS::FASTA& fasta,
        HDF5::File& hdf5,
        const MPI::Manager& mpi,
        const Params& params,
        Stats& stats
    );

    void run();

};

Mode mode(bool lowmem, bool joint);
Spread spread(bool uniform, bool mutation_informed);

std::unique_ptr<__Main> get_main(
    TinyHTS::File& file,
    TinyHTS::FASTA& fasta,
    HDF5::File& hdf5,
    const MPI::Manager& mpi,
    const Params& params,
    Stats& stats
);


} // namespace cmuts

#endif
