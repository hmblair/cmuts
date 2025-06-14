#include "cmuts.hpp"
#include "common.hpp"

namespace cmuts {





//
// Helpers for joint
//





template <typename dtype, bool match>
static inline void __joint(view_t<dtype, _ndims(Mode::Joint)> arr, int32_t ix, dtype mask) {

    if constexpr (match) {

        for (int32_t jx = 0; jx < ix; jx++) {

            dtype nval   = arr.periodic(jx, -1, 0, 0);
            dtype nval_c = 1 - nval;

            arr(ix, jx, NOMOD, MOD)   += mask * nval;
            arr(jx, ix, NOMOD, MOD)   += mask * nval;

            arr(ix, jx, NOMOD, NOMOD) += mask * nval_c;
            arr(jx, ix, NOMOD, NOMOD) += mask * nval_c;

        }

        arr(ix, ix, NOMOD, NOMOD) += mask;

    } else {

        for (int32_t jx = 0; jx < ix; jx++) {

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





static inline int32_t _first_match(
    base_t base,
    const seq_t& seq,
    int32_t M,
    int32_t N
) {

    for (int32_t i = M; i <= N; i++) {
        if (seq[i] == base) { return i+1; }
    }

    return -1;

}


static inline int32_t _get_ambiguous_end(
    int32_t start,
    int32_t end,
    const seq_t& sequence
) {

    auto size = static_cast<int32_t>(sequence.size());
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

    int32_t M = start;
    int32_t N = end;

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
    int32_t ix
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
        return arr(ix, LOWMEM_MOD);
    }

    if constexpr (mode == Mode::Joint) {
        return arr.periodic(ix, -1, 0, 0);
    }

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

template <typename dtype, Mode mode>
static inline std::vector<dtype> _spread_weights(
    int32_t start,
    int32_t end,
    view_t<dtype, _ndims(mode)> arr,
    dtype mask,
    Spread spread
) {

    int32_t length = end - start;
    std::vector<dtype> weights(length, 0);

    switch (spread) {

        case Spread::None: {

            weights[length - 1] = mask;
            break;

        }

        case Spread::Uniform: {

            for (int32_t ix = start; ix < end; ix++) {
                weights[ix - start] = mask / length;
            }

        }

        case Spread::MutationInformed: {

            for (int32_t ix = start; ix < end; ix++) {
                weights[ix - start] = _total_muts<dtype, mode>(arr, ix) * mask;
            }

            _normalize<dtype>(weights);
            break;

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
    int32_t rpos,
    base_t rbase,
    dtype mask,
    const Params& params
) {

    dtype _val = mask || !params.filter_coverage;

    // Count the base type and position
    if constexpr (mode == Mode::Normal) {
        arr(rpos, rbase, rbase) += _val;
    }

    // Count the base position
    if constexpr (mode == Mode::LowMem) {
        arr(rpos, LOWMEM_COV) += _val;
    }

    // Invoke the joint counter
    if constexpr (mode == Mode::Joint) {
        __joint<dtype, true>(arr, rpos, mask);
    }

};


template <typename dtype, Mode mode>
static inline void __mismatch_core(
    view_t<dtype, _ndims(mode)> arr,
    int32_t rpos,
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
        arr(rpos, LOWMEM_MOD) += mask;
        arr(rpos, LOWMEM_COV) += mask;
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
    int32_t rpos,
    base_t qbase,
    dtype mask
) {

    // Record the new base and position
    if constexpr (mode == Mode::Normal) {
        arr(rpos, qbase, IX_INS) += mask;
    }

    // Record the position; no associated match
    if constexpr (mode == Mode::LowMem) {
        arr(rpos, LOWMEM_MOD) += mask;
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
    int32_t rpos,
    base_t rbase,
    dtype mask
) {

    // Record the original base and position
    if constexpr (mode == Mode::Normal) {
        arr(rpos, rbase, IX_DEL) += mask;
    }

    // Record the position and an associated match
    if constexpr (mode == Mode::LowMem) {
        arr(rpos, LOWMEM_MOD) += mask;
        arr(rpos, LOWMEM_COV) += mask;
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
    HTS::CIGAR_op& op,
    int32_t& rpos,
    int32_t& qpos,
    const seq_t& reference,
    const std::vector<dtype>& mask,
    const Params& params
) {

    while (op.advance()) {

        if constexpr (CONSUMES_RPOS) { rpos--; }
        if constexpr (CONSUMES_QPOS) { qpos--; }
        base_t rbase = reference[rpos];

        __match_core<dtype, mode>(arr, rpos, rbase, mask[qpos], params);

    }

}


template <typename dtype, Mode mode>
static inline void __mismatch(
    view_t<dtype, _ndims(mode)> arr,
    HTS::CIGAR_op& op,
    int32_t& rpos,
    int32_t& qpos,
    int32_t& last,
    const seq_t& reference,
    const std::vector<dtype>& mask,
    const Params& params
) {

    rpos--;
    qpos--;
    op.advance();
    base_t rbase = reference[rpos];
    base_t qbase = op.last(rbase);

    if (params.mismatches && last - rpos >= params.collapse && qbase != IX_UNK) {
        __mismatch_core<dtype, mode>(arr, rpos, qbase, rbase, mask[qpos]);
        last = rpos;
    } else {
        __match_core<dtype, mode>(arr, rpos, rbase, mask[qpos], params);
    }

    __match<dtype, mode, true, true>(arr, op, rpos, qpos, reference, mask, params);

}


template <typename dtype, Mode mode>
static inline void __ins(
    view_t<dtype, _ndims(mode)> arr,
    HTS::CIGAR_op& op,
    int32_t& rpos,
    int32_t& qpos,
    int32_t& last,
    const seq_t& reference,
    const std::vector<dtype>& mask,
    const Params& params
) {

    if (rpos >= reference.size()) {
        qpos -= op.length();
        return;
    }

    qpos--;
    base_t qbase = op.last();
    if (params.insertions && last - rpos >= params.collapse && qbase != IX_UNK) {
        __ins_core<dtype, mode>(arr, rpos, qbase, mask[qpos]);
        last = rpos;
    }

    qpos -= (op.length() - 1);

}


template <typename dtype, Mode mode>
static inline void __del(
    view_t<dtype, _ndims(mode)> arr,
    HTS::CIGAR_op& op,
    int32_t& rpos,
    int32_t& qpos,
    int32_t& last,
    const seq_t& reference,
    const std::vector<dtype>& mask,
    const Params& params
) {

    rpos--;
    op.advance();

    int32_t start     = rpos - op.length() + 1;
    int32_t end       = rpos;
    int32_t ambig_end = _get_ambiguous_end(start, end, reference);

    // TODO: remove this and fix
    ambig_end = end + 1;

    std::vector<dtype> weights = _spread_weights<dtype, mode>(end, ambig_end, arr, mask[qpos], params.spread);

    for (int32_t ix = ambig_end - 1; ix >= end; ix--) {

        base_t rbase = reference[ix];
        if (params.deletions && last - ix >= params.collapse) {
            __del_core<dtype, mode>(arr, ix, rbase, weights[ix - end]);
            last = ix;
        } else {
            __match_core<dtype, mode>(arr, ix, rbase, weights[ix - end], params);
        }

    }

    __match<dtype, mode, true, false>(arr, op, rpos, qpos, reference, mask, params);

}


template <typename dtype, Mode mode>
static inline void __count(
    HTS::Alignment& aln,
    const seq_t& reference,
    view_t<dtype, _ndims(mode)> arr,
    const Params& params
) {

    HTS::CIGAR& cigar = aln.cigar;
    const HTS::PHRED& phred = aln.phred;

    // Initial reference and query positions

    int32_t rpos = aln.offset + cigar.rlength();
    int32_t qpos = aln.length;

    // Position of the most recent (3'-wise) mutation

    int32_t last = rpos + params.collapse;

    // PHRED quality mask

    std::vector<dtype> mask = phred.mask<dtype>(params.min_phred, params.quality_window);

    // Loop backwards (3' -> 5') along the CIGAR

    for (auto& op : std::ranges::reverse_view(cigar)) {

        switch (op.type()) {

            case HTS::CIGAR_t::MATCH: {

                __match<dtype, mode, true, true>(arr, op, rpos, qpos, reference, mask, params);
                break;

            }

            case HTS::CIGAR_t::MISMATCH: {

                __mismatch<dtype, mode>(arr, op, rpos, qpos, last, reference, mask, params);
                break;

            }

            case HTS::CIGAR_t::DEL: {

                if (op.length() > params.max_indel_length) {
                    __match<dtype, mode, true, false>(arr, op, rpos, qpos, reference, mask, params);
                    break;
                }

                __del<dtype, mode>(arr, op, rpos, qpos, last, reference, mask, params);
                break;

            }

            case HTS::CIGAR_t::INS: {

                if (op.length() > params.max_indel_length) {
                    qpos -= op.length();
                    break;
                }

                __ins<dtype, mode>(arr, op, rpos, qpos, last, reference, mask, params);
                break;

            }

            case HTS::CIGAR_t::SOFT: {

                qpos -= op.length();
                break;

            }

            case HTS::CIGAR_t::HARD: {

                break;

            }

            case HTS::CIGAR_t::SKIP: {

                rpos -= op.length();
                break;

            }

            case HTS::CIGAR_t::PAD:
            case HTS::CIGAR_t::BACK: {

                break;

            }

            case HTS::CIGAR_t::UNKNOWN: {

                throw std::runtime_error("Unknown CIGAR operation encountered at reference position " + std::to_string(rpos) + ".");

            }

        }

    }

}


static inline bool __check_quality(
    const HTS::Alignment& aln,
    const Params& params
) {
    return (
        aln.aligned                       &&
        aln.mapq    >= params.min_mapq    &&
        aln.length  >= params.min_length  &&
        aln.length  <= params.max_length
    );
}


std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<> dis(0.0, 1.0);
static inline float __sample() { return dis(gen); }

template <typename dtype, Mode mode, bool subsample>
static inline void __count_with_quality_check(
    HTS::Alignment& aln,
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
        __count<dtype, mode>(aln, reference, arr, params);
        stats.processed();
    }

}

template <typename dtype, Mode mode, bool subsample>
static inline void __count_reference(
    HTS::Iterator& iter,
    const seq_t& reference,
    view_t<dtype, _ndims(mode)> arr,
    const Params& params,
    Stats& stats
) {

    while (!iter.end()) {

        HTS::Alignment aln = iter.next();
        __count_with_quality_check<dtype, mode, subsample>(aln, reference, arr, params, stats);

        // Zero out the temp portion of the array
        if constexpr (mode == Mode::Joint) {
            xt::view(arr, xt::all(), -1, 0).fill(0.0);
        }

    }

}





//
// Main
//





static inline std::vector<size_t> __dims(const BinaryFASTA& fasta, Mode mode) {

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


template<typename dtype, size_t N>
static inline HDF5::Memspace<dtype, N> __memspace(
    BinaryFASTA& fasta,
    HDF5::File& hdf5,
    const std::string& name,
    Mode mode
) {

    std::vector<size_t> dims = __dims(fasta, mode);
    return hdf5.memspace<dtype, N>(dims, name);

}


template <typename dtype, Mode mode>
static inline void __fill_at_offset(
    HTS::Iterator& iter,
    const seq_t& reference,
    HDF5::Memspace<dtype, _ndims(mode)>& memspace,
    const Params& params,
    Stats& stats,
    int32_t offset
) {

    view_t<dtype, _ndims(mode)> arr = memspace.view(offset);
    if (params.subsample < 1.0) {
        return __count_reference<dtype, mode, true>(iter, reference, arr, params, stats);
    } else {
        return __count_reference<dtype, mode, false>(iter, reference, arr, params, stats);
    }

}


template <typename dtype, Mode mode>
static inline void __process(
    HTS::File& file,
    BinaryFASTA& fasta,
    HDF5::File& hdf5,
    const Params& params,
    Stats& stats,
    HDF5::Memspace<dtype, _ndims(mode)>& memspace,
    int32_t min
) {

    int32_t max = std::min(
        min + memspace.size(),
        file.size()
    );

    for (int32_t ix = min; ix < max; ix++) {

        bool seek = (ix % hdf5.chunk_size() == 0);
        std::shared_ptr<HTS::Iterator> iter = file.get(ix, seek);
        seq_t reference = fasta.sequence(ix);

        int32_t offset = ix - min;
        __fill_at_offset<dtype, mode>(*iter, reference, memspace, params, stats, offset);

    }

    memspace.safe_write(min);
    memspace.clear();

}


Main::Main(
    HTS::File& file,
    BinaryFASTA& fasta,
    HDF5::File& hdf5,
    const MPI::Manager& mpi,
    const Params& params,
    Stats& stats
) : file(file), fasta(fasta), hdf5(hdf5), mpi(mpi), params(params), chunk(mpi.chunk(hdf5.chunk_size(), file.size())), stats(stats) {};

template <typename dtype, Mode mode>
TemplatedMain<dtype, mode>::TemplatedMain(
    HTS::File& file,
    BinaryFASTA& fasta,
    HDF5::File& hdf5,
    const MPI::Manager& mpi,
    const Params& params,
    Stats& stats
) : Main(file, fasta, hdf5, mpi, params, stats),
    memspace(__memspace<dtype, _ndims(mode)>(fasta, hdf5, _path(file.name()), mode)) {}

template <typename dtype, Mode mode>
void TemplatedMain<dtype, mode>::run() {

    stats.body();

    for (int32_t ix = chunk.low; ix < chunk.high; ix += chunk.step) {

        __process<dtype, mode>(
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

template class TemplatedMain<float, Mode::Normal>;
template class TemplatedMain<float, Mode::LowMem>;
template class TemplatedMain<float, Mode::Joint>;

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

std::unique_ptr<Main> _get_main(
    HTS::File& file,
    BinaryFASTA& fasta,
    HDF5::File& hdf5,
    const MPI::Manager& mpi,
    const Params& params,
    Stats& stats
) {

    switch (params.mode) {

        case Mode::Normal: {

                return std::make_unique<TemplatedMain<float, Mode::Normal>>(
                    file, fasta, hdf5, mpi, params, stats
                );

       }

        case Mode::LowMem: {

                return std::make_unique<TemplatedMain<float, Mode::LowMem>>(
                    file, fasta, hdf5, mpi, params, stats
                );

       }

        case Mode::Joint: {

                return std::make_unique<TemplatedMain<float, Mode::Joint>>(
                    file, fasta, hdf5, mpi, params, stats
                );

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





Stats::Stats(int64_t aligned, int64_t unaligned, int32_t references, const MPI::Manager& mpi) 
    : _aligned(aligned), _unaligned(unaligned), _references(references), _mpi(mpi) {

        if (mpi.root()) {
            _processed += unaligned;
            _skipped   += unaligned;
        }

    }


void Stats::processed() {
    _processed++;
}


void Stats::processed(int64_t n) {
    _processed += n;
}


void Stats::skipped() {
    _processed++;
    _skipped++;
}


void Stats::skipped(int64_t n) {
    _processed += n;
    _skipped += n;
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

    _mpi.up(3);

    if (_mpi.root()) {
        double processed = _percent(_processed, _aligned + _unaligned);
        double skipped   = _percent(_skipped, _aligned + _unaligned);
        _print_processed.print(processed);
        _print_skipped.print(skipped);
        _print_elapsed.print(_mpi.time_str());
    }

    _mpi.divide();
    _mpi.up();

}





} // namespace cmuts
