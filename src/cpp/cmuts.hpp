#ifndef _CMUTS_HEADER
#define _CMUTS_HEADER

#include <random>
#include <ranges>

#include "common.hpp"
#include "fasta.hpp"
#include "bam.hpp"
#include "cram.hpp"
#include "hdf5.hpp"
#include "mpi.hpp"
#include "utils.hpp"

const int32_t N_PAIRS    = 2;
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





std::vector<bool> _ignore_str_to_bool(const std::string& ignore);





//
// Enums
//





enum class Mode {

    // Count mutation types, locations, and coverage
    // Size: (N, L, 4, 7)
    Normal,
    // Compute the joint distribution of modifications
    // Size: (N, L, L, 2, 2)
    Pairwise,
    // Tokenize the reference sequences
    // Size: (N, L)
    Tokenize

};


constexpr size_t _ndims(Mode mode) {

    switch (mode) {
        case Mode::Normal:   { return 4; }
        case Mode::Pairwise: { return 5; }
        case Mode::Tokenize: { return 2; }
    }

}


enum class Spread {

    None,
    Uniform,
    MutationInformed

};





//
// Data
//


template <typename dtype, Mode mode>
class DataView {
private:

    HDF5::Memspace<dtype, _ndims(mode)> _memspace;
    int32_t _offset = 0;

public:

    // view_t<dtype, _ndims(mode)> arr;

    DataView<dtype, mode>(HDF5::Memspace<dtype, _ndims(mode)> memspace);
    view_t<dtype, _ndims(mode)> view();
    void update(int32_t offset);
    void write(int32_t offset);
    int32_t size() const;

};


template <typename dtype>
class Data {
public:

    Data<dtype>(
        BinaryFASTA& fasta,
        HDF5::File& hdf5,
        const std::string& name,
        bool pairwise
    );

    DataView<dtype, Mode::Normal> mods;
    std::optional<DataView<dtype, Mode::Pairwise>> pairs;

    // Mutation profile of the current read; only used during --pairwise

    std::vector<dtype> tmp;
    int32_t min = 0;

    void update(int32_t offset);
    void write(int32_t offset);
    void size() const;

};





//
// Stats
//





class Stats {
private:

    int64_t _processed  = 0;
    int64_t _skipped    = 0;

    int64_t _bases_skipped   = 0;
    int64_t _bases_processed = 0;

    int64_t _aligned    = 0;
    int64_t _unaligned  = 0;

    int32_t _references = 0;
    int32_t _length     = 0;

    int64_t _files      = 0;
    int64_t _curr_file  = 0;

    const MPI::Manager& _mpi;

    Utils::Line _print_files     = Utils::Line("File");
    Utils::Line _print_processed = Utils::Line("Reads processed", "%");
    Utils::Line _print_skipped   = Utils::Line("Low-quality reads", "%");
    Utils::Line _print_bases_skipped   = Utils::Line("Low-quality bases", "%");
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


    void update_bases(int64_t skipped, int64_t total);

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

    bool pairwise;
    Spread spread;
    int32_t min_mapq;
    int32_t min_phred;
    int32_t min_length;
    int32_t max_length;
    int32_t max_indel_length;
    int32_t quality_window;
    int32_t collapse;
    int32_t max_hamming;
    bool mismatches;
    bool insertions;
    bool deletions;
    bool forward;
    bool reverse;
    float subsample;
    bool no_filter_matches;
    bool no_filter_insertions;
    bool no_filter_deletions;
    std::vector<bool> bases;
    bool ambiguous;
    int32_t gap;
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





template <typename dtype>
void run(
    HTS::File& file,
    BinaryFASTA& fasta,
    HDF5::File& hdf5,
    const MPI::Manager& mpi,
    const Params& params,
    Stats& stats,
    const std::string& name
);





Mode mode(bool lowmem, bool joint);
Spread spread(bool uniform, bool mutation_informed);





} // namespace cmuts





#endif // _CMUTS_HEADER
