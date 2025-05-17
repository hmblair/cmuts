#include "cmuts.hpp"
#include <stdexcept>

namespace cmuts {





//
// Helpers for joint
//





template <typename dtype, bool match>
static inline void __joint(view_t<dtype, _ndims(Mode::Joint)> arr, int64_t ix, dtype mask) {

    if constexpr (match) {

        for (int64_t jx = 0; jx < ix; jx++) {

            dtype nval   = arr.periodic(jx, -1, 0, 0);
            dtype nval_c = 1 - nval;

            arr(ix, jx, NOMOD, MOD)   += mask * nval;
            arr(jx, ix, NOMOD, MOD)   += mask * nval;

            arr(ix, jx, NOMOD, NOMOD) += mask * nval_c;
            arr(jx, ix, NOMOD, NOMOD) += mask * nval_c;

        }

        arr(ix, ix, NOMOD, NOMOD) += mask;

    } else {

        for (int64_t jx = 0; jx < ix; jx++) {

            dtype nval   = arr.periodic(jx, -1, 0, 0);
            dtype nval_c = 1 - nval;

            arr(ix, jx, MOD, MOD)   += mask * nval;
            arr(jx, ix, MOD, MOD)   += mask * nval;

            arr(ix, jx, MOD, NOMOD) += mask * nval_c;
            arr(jx, ix, MOD, NOMOD) += mask * nval_c;

        }

        arr(ix, ix, MOD, MOD) += mask;

    }

}



//
// Helpers for ambiguous deletions
//



static inline int64_t _first_match(
    base_t base,
    const seq_t& seq,
    int64_t M,
    int64_t N
) {

    for (int64_t i = M; i <= N; i++) {

        if (seq[i] == base) { return i+1; }

    }

    return -1;

}

static inline int64_t _get_ambiguous_end(
    int64_t start,
    int64_t end,
    const seq_t& sequence
) {

    int64_t size = sequence.size();
    if (end >= size) {
        return size;
    }

    // We wish to find the final base for which the deletion could
    // have occured at. Given a string S, with deletion starting at
    // position M and ending at position N, we wish to find the largest
    // K such that S[K] could potentially be part of the deletion.
    // If D is the string with the deletion, then we make use of
    // the fact that S[K] is a valid deletion if and only if 
    // D[M:K] is a subsequence of S[M:K-1].

    // We can find the largest such K in an efficient manner.
    //    - Check if D[N] is in S[M:N]. While this is so,
    //    - Advance M to the base after the first position in
    //      S where D[N] occurs, and
    //    - Update N -> N + 1.
    // The final N is the furthest base in the 3' direction
    // that can be deleted to yield the same deletion. If N
    // ever reaches the end of the reference sequence, we
    // must terminate.

    int64_t M = start, N = end;

    while (M != -1 && N != size) {
        M = _first_match(sequence[N + 1], sequence, M, N);
        N++;
    }

    // This is the position of the base _after_ the final base that can be deleted.
    return N;

}



//
// Helpers for spreading deletions
//



template <typename dtype, Mode mode>
static inline dtype _total_muts(
    view_t<dtype, _ndims(mode)> arr,
    int64_t ix
) {

    if constexpr (mode == Mode::Normal) {
        dtype total = 0;
        // From A
        total += arr(ix, IX_A, IX_C);
        total += arr(ix, IX_A, IX_G);
        total += arr(ix, IX_A, IX_T);
        // From C
        total += arr(ix, IX_C, IX_A);
        total += arr(ix, IX_C, IX_G);
        total += arr(ix, IX_C, IX_T);
        // From G
        total += arr(ix, IX_G, IX_A);
        total += arr(ix, IX_G, IX_C);
        total += arr(ix, IX_G, IX_T);
        // From T
        total += arr(ix, IX_T, IX_A);
        total += arr(ix, IX_T, IX_C);
        total += arr(ix, IX_T, IX_G);
        return total;
    }

    if constexpr (mode == Mode::LowMem) {
        return arr(ix, 0);
    }

    if constexpr (mode == Mode::Joint) {
        return arr.periodic(ix, -1, 0, 0);
    }

    return static_cast<dtype>(0);

}

template <typename dtype>
static inline void _normalize(std::vector<dtype>& arr) {

    dtype total = 0;
    for (const dtype& val : arr) { total += val; }

    if (total == 0) {

        for (dtype& val : arr) { val = 1 / arr.size(); }

    } else {

        for (dtype& val : arr) { val /= total; }

    }

}

template <typename dtype, Mode mode, Spread spread>
static inline std::vector<dtype> _spread_weights(
    int64_t start,
    int64_t end,
    view_t<dtype, _ndims(mode)> arr,
    dtype mask
) {

    int64_t length = end - start;
    std::vector<dtype> weights(length, 0);

    if constexpr (spread == Spread::None) {

        weights[length - 1] = mask;

    }

    if constexpr (spread == Spread::Uniform) {

        for (int64_t ix = start; ix < end; ix++) {
            weights[ix - start] = _total_muts<dtype, mode>(arr, ix) * mask;
        }

        _normalize<dtype>(weights);

    }

    if constexpr (spread == Spread::Uniform) {

        for (int64_t ix = start; ix < end; ix++) {
            weights[ix - start] = mask / length;
        }

    }

    return weights;

}



//
// Mutation counting
//



template <typename dtype, Mode mode>
static inline void __match_core(
    view_t<dtype, _ndims(mode)> arr,
    int64_t rpos,
    base_t rbase,
    dtype mask
) {

    // Count the base type and position
    if constexpr (mode == Mode::Normal) {
        arr(rpos, rbase, rbase) += mask;
    }

    // Count the base position
    if constexpr (mode == Mode::LowMem) {
        arr(rpos, 1) += mask;
    }

    // Invoke the joint counter
    if constexpr (mode == Mode::Joint) {
        __joint<dtype, true>(arr, rpos, mask);
    }

};

template <typename dtype, Mode mode>
static inline void __mismatch_core(
    view_t<dtype, _ndims(mode)> arr,
    int64_t rpos,
    base_t qbase,
    base_t rbase,
    dtype mask
) {

    // Record the original base, new base, and position
    if constexpr (mode == Mode::Normal) {
        arr(rpos, rbase, qbase) += mask;
    }

    // Record the position and an associated match
    if constexpr (mode == Mode::LowMem) {
        arr(rpos, 0) += mask;
        arr(rpos, 1) += mask;
    }

    // Store the modification in the temp dim and invoke the joint counter
    if constexpr (mode == Mode::Joint) {
        arr.periodic(rpos, -1, 0, 0) += mask;
        __joint<dtype, false>(arr, rpos, mask);
    }

};

template <typename dtype, Mode mode>
static inline void __ins_core(
    view_t<dtype, _ndims(mode)> arr,
    int64_t rpos,
    base_t qbase,
    dtype mask
) {

    // Record the new base and position
    if constexpr (mode == Mode::Normal) {
        arr(rpos, qbase, IX_INS) += mask;
    }

    // Record the position; no associated match
    if constexpr (mode == Mode::LowMem) {
        arr(rpos, 0) += mask;
    }

    // Store the modification in the temp dim and invoke the joint counter
    if constexpr (mode == Mode::Joint) {
        arr.periodic(rpos, -1, 0) += mask;
        __joint<dtype, false>(arr, rpos, mask);
    }

};

template <typename dtype, Mode mode>
static inline void __del_core(
    view_t<dtype, _ndims(mode)> arr,
    int64_t rpos,
    base_t rbase,
    dtype mask
) {

    // Record the original base and position
    if constexpr (mode == Mode::Normal) {
        arr(rpos, rbase, IX_DEL) += mask;
    }

    // Record the position and an associated match
    if constexpr (mode == Mode::LowMem) {
        arr(rpos, 0) += mask;
        arr(rpos, 1) += mask;
    }

    // Store the modification in the temp dim and invoke the joint counter
    if constexpr (mode == Mode::Joint) {
        arr.periodic(rpos, -1, 0, 0) += mask;
        __joint<dtype, false>(arr, rpos, mask);
    }

};

template <
    typename dtype,
    Mode mode,
    bool CONSUMES_RPOS,
    bool CONSUMES_QPOS
>
static inline void __match(
    view_t<dtype, _ndims(mode)> arr,
    TinyHTS::CIGAR_op& op,
    int64_t& rpos,
    int64_t& qpos,
    const TinyHTS::Alignment& aln,
    const seq_t& reference,
    const std::vector<dtype>& mask
) {

    while (op.advance()) {

        if constexpr (CONSUMES_RPOS) { rpos--; }
        if constexpr (CONSUMES_QPOS) { qpos--; }
        base_t rbase = reference[rpos];

        __match_core<dtype, mode>(arr, rpos, rbase, mask[qpos]);

    }

}

template <typename dtype, Mode mode>
static inline void __mismatch(
    view_t<dtype, _ndims(mode)> arr,
    TinyHTS::CIGAR_op& op,
    int64_t& rpos,
    int64_t& qpos,
    int64_t& last,
    const TinyHTS::Alignment& aln,
    const seq_t& reference,
    const std::vector<dtype>& mask,
    const Params& params
) {

    rpos--;
    qpos--;
    op.advance();
    base_t rbase = reference[rpos];
    base_t qbase = aln.base(qpos);

    if (params.mismatches && last - rpos >= params.collapse && qbase != IX_UNK) {
        __mismatch_core<dtype, mode>(arr, rpos, qbase, rbase, mask[qpos]);
        last = rpos;
    } else {
        __match_core<dtype, mode>(arr, rpos, rbase, mask[qpos]);
    }

    __match<dtype, mode, true, true>(arr, op, rpos, qpos, aln, reference, mask);

}

template <typename dtype, Mode mode>
static inline void __ins(
    view_t<dtype, _ndims(mode)> arr,
    TinyHTS::CIGAR_op& op,
    int64_t& rpos,
    int64_t& qpos,
    int64_t& last,
    const TinyHTS::Alignment& aln,
    const seq_t& reference,
    const std::vector<dtype>& mask,
    const Params& params
) {

    if (rpos >= reference.size()) {
        qpos -= op.length();
        return;
    }

    qpos--;
    base_t qbase = aln.base(qpos);
    if (params.insertions && last - rpos >= params.collapse && qbase != IX_UNK) {
        __ins_core<dtype, mode>(arr, rpos, qbase, mask[qpos]);
        last = rpos;
    }

    qpos -= (op.length() - 1);

}

template <typename dtype, Mode mode, Spread spread>
static inline void __del(
    view_t<dtype, _ndims(mode)> arr,
    TinyHTS::CIGAR_op& op,
    int64_t& rpos,
    int64_t& qpos,
    int64_t& last,
    const TinyHTS::Alignment& aln,
    const seq_t& reference,
    const std::vector<dtype>& mask,
    const Params& params
) {

    rpos--;
    op.advance();

    int64_t start     = rpos - op.length() + 1;
    int64_t end       = rpos;
    int64_t ambig_end = _get_ambiguous_end(start, end, reference);

    // TODO: remove this and fix
    ambig_end = end + 1;

    std::vector<dtype> weights = _spread_weights<dtype, mode, spread>(end, ambig_end, arr, mask[qpos]);

    for (int64_t ix = ambig_end - 1; ix >= end; ix--) {

        base_t rbase = reference[ix];
        if (params.deletions && last - ix >= params.collapse) {
            __del_core<dtype, mode>(arr, ix, rbase, weights[ix - end]);
            last = ix;
        } else {
            __match_core<dtype, mode>(arr, ix, rbase, weights[ix - end]);
        }

    }

    __match<dtype, mode, true, false>(arr, op, rpos, qpos, aln, reference, mask);

}

template <typename dtype, Mode mode, Spread spread>
static inline void __count(
    const TinyHTS::Alignment& aln,
    const seq_t& reference,
    view_t<dtype, _ndims(mode)> arr,
    const Params& params
) {

    // Matches, mismatches, etc.
    TinyHTS::CIGAR cigar = aln.cigar();
    // Initial reference and query positions
    int64_t rpos = aln.offset() + cigar.rlength(), qpos = aln.length();
    // Position of the most recent (3'-wise) mutation
    int64_t last = rpos + params.collapse;
    // PHRED quality mask
    std::vector<dtype> mask = aln.mask<dtype>(params.min_phred, params.quality_window);

    // Loop backwards (3' -> 5') along the CIGAR

    for (auto& op : std::ranges::reverse_view(cigar)) {

        switch (op.type()) {

            case TinyHTS::CIGAR_t::MATCH: {

                __match<dtype, mode, true, true>(arr, op, rpos, qpos, aln, reference, mask);
                break;

            }

            case TinyHTS::CIGAR_t::MISMATCH: {

                __mismatch<dtype, mode>(arr, op, rpos, qpos, last, aln, reference, mask, params);
                break;

            }

            case TinyHTS::CIGAR_t::DEL: {

                if (op.length() > params.max_indel_length) {
                    __match<dtype, mode, true, false>(arr, op, rpos, qpos, aln, reference, mask);
                    break;
                }

                __del<dtype, mode, spread>(arr, op, rpos, qpos, last, aln, reference, mask, params);
                break;

            }

            case TinyHTS::CIGAR_t::INS: {

                if (op.length() > params.max_indel_length) {
                    qpos -= op.length();
                    break;
                }

                __ins<dtype, mode>(arr, op, rpos, qpos, last, aln, reference, mask, params);
                break;

            }

            case TinyHTS::CIGAR_t::SOFT: {

                qpos -= op.length();
                break;

            }

            case TinyHTS::CIGAR_t::HARD: {

                break;

            }

            case TinyHTS::CIGAR_t::SKIP: {

                rpos -= op.length();
                break;

            }

            case TinyHTS::CIGAR_t::PAD: {

                break;

            }

            case TinyHTS::CIGAR_t::BACK: {

                break;

            }

            case TinyHTS::CIGAR_t::UNKNOWN: {

                throw std::runtime_error("Unknown CIGAR operation encountered at reference position " + std::to_string(rpos) + ".");

            }

        }

    }

}

static inline bool __check_quality(
    const TinyHTS::Alignment& aln,
    const Params& params
) {
    return (
        aln.aligned() &&
        aln.mapq()    >= params.min_mapq   &&
        aln.length() >= params.min_length &&
        aln.length() <=  params.max_length
    );
}

std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<> dis(0.0, 1.0);
static inline float __sample() { return dis(gen); }

template <typename dtype, Mode mode, Spread spread, bool subsample>
static inline void __count_with_quality_check(
    const TinyHTS::Alignment& aln,
    const seq_t& reference,
    view_t<dtype, _ndims(mode)> arr,
    const Params& params,
    Stats& stats
) {

    if constexpr (subsample) {
        if (params.subsample < __sample()) { return; }
    }

    if (!__check_quality(aln, params)) {
        stats.skipped();
    } else {
        __count<dtype, mode, spread>(aln, reference, arr, params);
        stats.processed();
    }

}

template <typename dtype, Mode mode, Spread spread, bool subsample>
static inline void __count_reference(
    TinyHTS::Alignment& aln,
    const seq_t& reference,
    view_t<dtype, _ndims(mode)> arr,
    const Params& params,
    Stats& stats
) {

    while (aln.next()) {

        __count_with_quality_check<dtype, mode, spread, subsample>(aln, reference, arr, params, stats);

        // Zero out the temp portion of the array
        if constexpr (mode == Mode::Joint) {
            xt::view(arr, xt::all(), -1, 0).fill(0.0);
        }

    }

}



//
// Main
//



static inline std::vector<size_t> __dims(const TinyHTS::FASTA& fasta, Mode mode) {

    size_t size   = fasta.size();
    size_t length = fasta.longest();

    switch (mode) {

        case Mode::Normal: {
            return {size, length, N_BASES, N_DELBASES};
        }

        case Mode::LowMem: {
            return {size, length, 2};
        }

        case Mode::Joint: {
            return {size, length, length + 1, 2, 2};
        }

        case Mode::Tokenize: {
            return {size, length};
        }

    }

}


static inline std::string _path(const std::string& name) {

    std::filesystem::path path(name);
    std::string stem = path.stem().string();
    std::string parent = path.parent_path().string();
    if (parent.empty()) { return stem; }
    return parent + "/" + stem;

}


template<typename dtype, size_t N>
static inline HDF5::Memspace<dtype, N> __memspace(
    TinyHTS::FASTA& fasta,
    HDF5::File& hdf5,
    const std::string& name,
    Mode mode
) {

    std::vector<size_t> dims = __dims(fasta, mode);
    return hdf5.memspace<dtype, N>(dims, name);

}

template <typename dtype, Mode mode, Spread spread>
static inline void __fill_at_offset(
    TinyHTS::Alignment& aln,
    const seq_t& reference,
    HDF5::Memspace<dtype, _ndims(mode)>& memspace,
    const Params& params,
    Stats& stats,
    int64_t offset
) {

    view_t<dtype, _ndims(mode)> arr = memspace.view(offset);
    if (params.subsample < 1.0) {
        return __count_reference<dtype, mode, spread, true>(aln, reference, arr, params, stats);
    } else {
        return __count_reference<dtype, mode, spread, false>(aln, reference, arr, params, stats);
    }

}

template <typename dtype, Mode mode, Spread spread>
static inline void __process(
    TinyHTS::File& file,
    TinyHTS::FASTA& fasta,
    HDF5::File& hdf5,
    const Params& params,
    Stats& stats,
    HDF5::Memspace<dtype, _ndims(mode)>& memspace,
    int64_t min
) {

    int64_t max = std::min(
        min + memspace.size(),
        file.size()
    );

    for (int64_t ix = min; ix < max; ix++) {

        bool seek = (ix % hdf5.chunk_size() == 0);
        TinyHTS::Alignment aln = file.alignment(ix, seek);
        seq_t reference = fasta.sequence(ix);

        int64_t offset = ix - min;
        __fill_at_offset<dtype, mode, spread>(aln, reference, memspace, params, stats, offset);

    }

    memspace.safe_write(min);
    memspace.clear();

}

__Main::__Main(
    TinyHTS::File& file,
    TinyHTS::FASTA& fasta,
    HDF5::File& hdf5,
    const MPI::Manager& mpi,
    const Params& params,
    Stats& stats
) : file(file), fasta(fasta), hdf5(hdf5), mpi(mpi), params(params), chunk(mpi.chunk(hdf5.chunk_size(), file.size())), stats(stats) { };

template <typename dtype, Mode mode, Spread spread>
Main<dtype, mode, spread>::Main(
    TinyHTS::File& file,
    TinyHTS::FASTA& fasta,
    HDF5::File& hdf5,
    const MPI::Manager& mpi,
    const Params& params,
    Stats& stats
) : __Main(file, fasta, hdf5, mpi, params, stats), memspace(__memspace<dtype, _ndims(mode)>(fasta, hdf5, _path(file.name()), mode)) {}

template <typename dtype, Mode mode, Spread spread>
void Main<dtype, mode, spread>::run() {

    stats.body();

    for (int64_t ix = chunk.low; ix < chunk.high; ix += chunk.step) {

        __process<dtype, mode, spread>(
            file,
            fasta,
            hdf5,
            params,
            stats,
            memspace,
            ix
        );
        stats.aggregate();
        stats.body();

    }

}

template class Main<float, Mode::Normal, Spread::None>;
template class Main<float, Mode::Normal, Spread::Uniform>;
template class Main<float, Mode::Normal, Spread::MutationInformed>;

template class Main<float, Mode::LowMem, Spread::None>;
template class Main<float, Mode::LowMem, Spread::Uniform>;
template class Main<float, Mode::LowMem, Spread::MutationInformed>;

template class Main<float, Mode::Joint, Spread::None>;
template class Main<float, Mode::Joint, Spread::Uniform>;
template class Main<float, Mode::Joint, Spread::MutationInformed>;

Mode mode(bool lowmem, bool joint) {

    if (lowmem && joint) {
        throw std::runtime_error("--low-mem and --joint are mutully exclusive.");
    }

    if (lowmem) { return Mode::LowMem; }
    if (joint) { return Mode::Joint; }
    return Mode::Normal;

}

Spread spread(bool uniform, bool mutation_informed) {

    if (uniform && mutation_informed) {
        throw std::runtime_error("--uniform-spread and --mutation-spread are mutually exclusive.");
    }

    if (uniform) { return Spread::Uniform; }
    if (mutation_informed) { return Spread::MutationInformed; }
    return Spread::None;

}

std::unique_ptr<__Main> get_main(
    TinyHTS::File& file,
    TinyHTS::FASTA& fasta,
    HDF5::File& hdf5,
    const MPI::Manager& mpi,
    const Params& params,
    Stats& stats
) {

    switch (params.mode) {

        case Mode::Normal: {
            switch (params.spread) {
                case Spread::None: {
                    return std::make_unique<Main<float, Mode::Normal, Spread::None>>(file, fasta, hdf5, mpi, params, stats);
                }
                case Spread::Uniform: {
                    return std::make_unique<Main<float, Mode::Normal, Spread::Uniform>>(file, fasta, hdf5, mpi, params, stats);
                }
                case Spread::MutationInformed: {
                    return std::make_unique<Main<float, Mode::Normal, Spread::MutationInformed>>(file, fasta, hdf5, mpi, params, stats);
                }
            }
        }

        case Mode::LowMem: {
            switch (params.spread) {
                case Spread::None: {
                    return std::make_unique<Main<float, Mode::LowMem, Spread::None>>(file, fasta, hdf5, mpi, params, stats);
                }
                case Spread::Uniform: {
                    return std::make_unique<Main<float, Mode::LowMem, Spread::Uniform>>(file, fasta, hdf5, mpi, params, stats);
                }
                case Spread::MutationInformed: {
                    return std::make_unique<Main<float, Mode::LowMem, Spread::MutationInformed>>(file, fasta, hdf5, mpi, params, stats);
                }
            }
        }

        case Mode::Joint: {
            switch (params.spread) {
                case Spread::None: {
                    return std::make_unique<Main<float, Mode::Joint, Spread::None>>(file, fasta, hdf5, mpi, params, stats);
                }
                case Spread::Uniform: {
                    return std::make_unique<Main<float, Mode::Joint, Spread::Uniform>>(file, fasta, hdf5, mpi, params, stats);
                }
                case Spread::MutationInformed: {
                    return std::make_unique<Main<float, Mode::Joint, Spread::MutationInformed>>(file, fasta, hdf5, mpi, params, stats);
                }
            }
        }

    }

}





//
// Stats
//





static inline double _percent(int64_t a, int64_t b) {
    double a_double = static_cast<double>(a);
    double b_double = static_cast<double>(b);
    return a_double / b_double * 100;
}

static inline void _print_header(int64_t aligned, int64_t unaligned, int64_t references) {

    Utils::Line _print_references("References");
    Utils::Line _print_aligned("Aligned reads");
    Utils::Line _print_unaligned("Unaligned reads");

    _print_references.print(references);
    _print_aligned.print(aligned);
    _print_unaligned.print(unaligned);
    Utils::divider();
    Utils::cursor_down(3);
    Utils::cursor_down(2);
    Utils::cursor_up(2);

}





//
// Stats
//





Stats::Stats(
    int64_t aligned,
    int64_t unaligned,
    int64_t references,
    const MPI::Manager& mpi
) : _aligned(aligned),
    _unaligned(unaligned),
    _references(references),
    _mpi(mpi) {}

void Stats::processed() {
    _processed++;
}

void Stats::skipped() {
    _processed++;
    _skipped++;
}

void Stats::aggregate() {
    _processed = _mpi.reduce(_processed);
    _skipped = _mpi.reduce(_skipped);
}

void Stats::header() const {
    if (_mpi.root()) {
        _print_header(_aligned, _unaligned, _references);
    }
}

void Stats::body() const {
    if (_mpi.root()) {
        Utils::cursor_up(3);
        double percent_processed = _percent(_processed, _aligned);
        double percent_skipped = _percent(_skipped, _aligned);
        _print_processed.print(percent_processed);
        _print_skipped.print(percent_skipped);
        _print_elapsed.print(_mpi.time_str());
    }
}



} // namespace cmuts
