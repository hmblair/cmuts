#ifndef _CMUTS_HEADER
#define _CMUTS_HEADER

#include "common.hpp"
#include "fasta.hpp"
#include "bam.hpp"
#include "cram.hpp"
#include "hdf5.hpp"
#include "mpi.hpp"
#include "utils.hpp"
#include <random>

const int32_t N_BASES    = 4;
const int32_t N_DELBASES = 7;

const std::string MODIFICATION_DS = "modifications";
const std::string JOINT_DS        = "joint";

const hts_pos_t NOMOD = 0;
const hts_pos_t MOD   = 1;

const hts_pos_t MOD_MOD     = 0;
const hts_pos_t MOD_NOMOD   = 1;
const hts_pos_t NOMOD_MOD   = 2;
const hts_pos_t NOMOD_NOMOD = 3;

const int32_t LOWMEM_MOD = 0;
const int32_t LOWMEM_COV = 1;

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
    int32_t _references = 0;
    int32_t _length     = 0;
    int64_t _files      = 0;
    int64_t _curr_file  = 0;

    const MPI::Manager& _mpi;

    Utils::Line _print_files     = Utils::Line("File");
    Utils::Line _print_processed = Utils::Line("Reads processed", "%");
    Utils::Line _print_skipped   = Utils::Line("Reads skipped", "%");
    Utils::Line _print_elapsed   = Utils::Line("Time elapsed");

public:

    Stats(
        int64_t files,
        int64_t aligned,
        int64_t unaligned,
        int32_t references,
        int32_t length,
        const MPI::Manager& mpi
    );

    void processed();
    void processed(int64_t n);
    void skipped();
    void skipped(int64_t n);
    void file();
    void aggregate();
    void header() const;
    void body() const;
    bool mod(int64_t n) const;

};





//
// Params
//





class Params {
public:

    Mode mode;
    Spread spread;
    int32_t min_mapq;
    int32_t min_phred;
    int32_t min_length;
    int32_t max_length;
    int32_t max_indel_length;
    int32_t quality_window;
    int32_t collapse;
    bool mismatches;
    bool insertions;
    bool deletions;
    float subsample;
    bool filter_coverage;
    bool ambiguous;
    bool contiguous;
    // std::vector<std::vector<bool>> valid;
    int64_t print_every;

};





//
// Ambiguous deletion detection
//





int32_t _get_ambiguous_end(
    int32_t start,
    int32_t end,
    const seq_t& sequence
);





//
// Main
//





class Main {
protected:

    HTS::File& file;
    BinaryFASTA& fasta;
    HDF5::File& hdf5;
    const MPI::Manager& mpi;
    const Params& params;
    MPI::Chunk chunk;
    Stats& stats;

public:

    Main(
        HTS::File& file,
        BinaryFASTA& fasta,
        HDF5::File& hdf5,
        const MPI::Manager& mpi,
        const Params& params,
        Stats& stats
    );
    virtual ~Main() = default;

    virtual void run() = 0;

};


template <typename dtype, Mode mode>
class TemplatedMain : public Main {
private:

    HDF5::Memspace<dtype, _ndims(mode)> memspace;

public:

    TemplatedMain(
        HTS::File& file,
        BinaryFASTA& fasta,
        HDF5::File& hdf5,
        const MPI::Manager& mpi,
        const Params& params,
        Stats& stats,
        const std::string& name
    );

    void run() override;

};


Mode mode(bool lowmem, bool joint);
Spread spread(bool uniform, bool mutation_informed);

std::unique_ptr<Main> _get_main(
    HTS::File& file,
    BinaryFASTA& fasta,
    HDF5::File& hdf5,
    const MPI::Manager& mpi,
    const Params& params,
    Stats& stats,
    const std::string& name
);





} // namespace cmuts





#endif // _CMUTS_HEADER
