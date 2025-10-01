#include "cmuts.hpp"
#include "common.hpp"

namespace cmuts {





//
// Helpers for joint
//



static int32_t _MAX_ = 0;
static std::vector<float> ARR_TMP;


template <typename dtype, bool match>
static inline void __joint(
    view_t<dtype, _ndims(Mode::Joint)> arr,
    std::vector<dtype> tmp,
    int32_t ix,
    dtype mask
) {

    if constexpr (match) {

        for (int32_t jx = _MAX_ - 1; jx > ix; jx--) {

            dtype nval   = tmp[jx];
            dtype nval_c = 1 - nval;

            arr(ix, jx, NOMOD, MOD)   += mask * nval;
            arr(jx, ix, MOD, NOMOD)   += mask * nval;

            arr(ix, jx, NOMOD, NOMOD) += mask * nval_c;
            arr(jx, ix, NOMOD, NOMOD) += mask * nval_c;

        }

        arr(ix, ix, NOMOD, NOMOD) += (mask * mask);
        arr(ix, ix, MOD, MOD)     += (1 - (mask * mask));

    } else {

        for (int32_t jx = _MAX_ - 1; jx > ix; jx--) {

            dtype nval   = tmp[jx];
            dtype nval_c = 1 - nval;

            arr(ix, jx, MOD, MOD)   += mask * nval;
            arr(jx, ix, MOD, MOD)   += mask * nval;

            arr(ix, jx, MOD, NOMOD) += mask * nval_c;
            arr(jx, ix, NOMOD, MOD) += mask * nval_c;

        }

        arr(ix, ix, MOD, MOD)     += (mask * mask);
        arr(ix, ix, NOMOD, NOMOD) += (1 - (mask * mask));

    }

}





//
// For ambiguous deletions
//





int32_t _get_ambiguous_end_contiguous(
    int32_t start,
    int32_t end,
    const seq_t& sequence
) {

    auto size = static_cast<int32_t>(sequence.size());
    if (end >= size) { return size; }

    int32_t M = start;
    int32_t N = end;

    while (N < size && sequence[M] == sequence[N]) { M++; N++; }

    return N;

}





int32_t _get_ambiguous_end(
    int32_t start,
    int32_t end,
    const seq_t& sequence
) {

    auto size = static_cast<int32_t>(sequence.size());
    if (end >= size) { return size; }

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

    while (M < N && N < size) {
        if (sequence[M] == sequence[N]) { N++; }
        M++;
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

    throw std::runtime_error("Invalid mode for _total_muts.");

}

template <typename dtype>
static inline void _normalize(std::vector<dtype>& arr, dtype mask) {

    dtype total = 0;
    for (const dtype& val : arr) { total += val; }

    if (total == 0) {
        for (dtype& val : arr) { val = mask / arr.size(); }

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
    if (length == 0) { return weights; }

    switch (spread) {

        case Spread::None: {

            weights[length - 1] = mask;
            break;

        }

        case Spread::Uniform: {

            for (int32_t ix = start; ix < end; ix++) {
                weights[ix - start] = mask / length;
            }
            break;

        }

        case Spread::MutationInformed: {

            for (int32_t ix = start; ix < end; ix++) {
                weights[ix - start] = _total_muts<dtype, mode>(arr, ix) * mask;
            }

            _normalize<dtype>(weights, mask);
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

    dtype _val = mask || params.no_filter_matches;

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
        ARR_TMP[rpos] = 0;
        __joint<dtype, true>(arr, ARR_TMP, rpos, _val);
    }

};


template <typename dtype, Mode mode>
static inline void __mismatch_core(
    view_t<dtype, _ndims(mode)> arr,
    int32_t rpos,
    base_t qbase,
    base_t rbase,
    dtype mask,
    const Params& params
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
        ARR_TMP[rpos] = mask;
        __joint<dtype, false>(arr, ARR_TMP, rpos, mask);
    }

};

template <typename dtype, Mode mode>
static inline void __ins_core(
    view_t<dtype, _ndims(mode)> arr,
    int32_t rpos,
    base_t qbase,
    dtype mask,
    const Params& params
) {

    dtype _val = mask || params.no_filter_insertions;

    // Record the new base and position
    if constexpr (mode == Mode::Normal) {
        arr(rpos, qbase, IX_INS) += _val;
    }

    // Record the position; no associated match
    if constexpr (mode == Mode::LowMem) {
        arr(rpos, LOWMEM_MOD) += _val;
    }

    // Store the modification in the temp dim and invoke the joint counter
    if constexpr (mode == Mode::Joint) {
        ARR_TMP[rpos] = mask;
        __joint<dtype, false>(arr, ARR_TMP, rpos, _val);
    }

};

template <typename dtype, Mode mode>
static inline void __del_core(
    view_t<dtype, _ndims(mode)> arr,
    int32_t rpos,
    base_t rbase,
    dtype mask,
    const Params& params
) {

    dtype _val = mask;

    // Record the original base and position
    if constexpr (mode == Mode::Normal) {
        arr(rpos, rbase, IX_DEL) += _val;
    }

    // Record the position and an associated match
    if constexpr (mode == Mode::LowMem) {
        arr(rpos, LOWMEM_MOD) += _val;
        arr(rpos, LOWMEM_COV) += _val;
    }

    // Store the modification in the temp dim and invoke the joint counter
    if constexpr (mode == Mode::Joint) {
        ARR_TMP[rpos] = mask;
        __joint<dtype, false>(arr, ARR_TMP, rpos, _val);
    }

};


template <typename dtype, Mode mode>
static inline void __term_core(
    view_t<dtype, _ndims(mode)> arr,
    int32_t rpos,
    base_t rbase,
    dtype mask
) {

    // Record the original base and position
    if constexpr (mode == Mode::Normal) {
        arr(rpos, rbase, IX_TERM) += mask;
    }

    // Store the modification in the temp dim and invoke the joint counter
    if constexpr (mode == Mode::Joint) {
        // arr.periodic(rpos, -1, 0, 0) = mask;
        // __joint<dtype, false>(arr, rpos, mask);
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
        __mismatch_core<dtype, mode>(arr, rpos, qbase, rbase, mask[qpos], params);
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

    // Insertions at the beginning of the sequence are ignored

    if (rpos > reference.size() || rpos <= 0) {
        qpos -= op.length();
        return;
    }

    qpos--;
    base_t qbase = op.last();
    if (params.insertions && last - rpos >= params.collapse && qbase != IX_UNK) {
        __ins_core<dtype, mode>(arr, rpos - 1, qbase, mask[qpos], params);
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

    // Skip conditions

    if (
        !params.deletions                     ||
        op.length() > params.max_indel_length ||
        last - (rpos - 1) < params.collapse
    ) {
        __match<dtype, mode, true, false>(arr, op, rpos, qpos, reference, mask, params);
        return;
    }

    // Start and end indices of the deletion

    int32_t start = rpos - op.length();
    int32_t end   = rpos;

    // Start and end indices of the ambiguous region
    // Defaults to simply the 3' base, but may change if ambiguous detection
    // is not disabled

    int32_t ambig_start = end - 1;
    int32_t ambig_end   = end;

    if (params.ambiguous) {
        if (params.contiguous) {
            ambig_end = _get_ambiguous_end_contiguous(start, end, reference);
        } else {
            ambig_end = _get_ambiguous_end(start, end, reference);
        }
    }

    dtype _val = mask[qpos] || params.no_filter_deletions;
    std::vector<dtype> weights = _spread_weights<dtype, mode>(ambig_start, ambig_end, arr, _val, params.spread);

    for (int32_t ix = ambig_end - 1; ix >= ambig_start; ix--) {
        base_t rbase = reference[ix];
        __del_core<dtype, mode>(arr, ix, rbase, weights[ix - ambig_start], params);
    }

    op.advance();
    rpos--;
    last = rpos;

    __match<dtype, mode, true, false>(arr, op, rpos, qpos, reference, mask, params);

}


template <typename dtype, Mode mode>
static inline void __count(
    HTS::Alignment& aln,
    const seq_t& reference,
    view_t<dtype, _ndims(mode)> arr,
    const Params& params
) {

    HTS::CIGAR& cigar       = aln.cigar;
    const HTS::PHRED& phred = aln.phred;

    // Initial reference and query positions

    int32_t rpos = aln.offset + cigar.rlength();
    int32_t qpos = aln.length;

    if constexpr(mode == Mode::Joint) {
        _MAX_ = rpos;
        ARR_TMP = std::vector<dtype>(reference.size(), 0);
    }

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

            case HTS::CIGAR_t::TERM: {

                __term_core<dtype, mode>(arr, rpos, reference[rpos], mask[qpos]);
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

    // Count termination events

    __term_core<dtype, mode>(arr, rpos, reference[rpos], mask[qpos]);

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

    HTS::Alignment aln;

    while (!iter.end()) {

        aln = iter.next();
        __count_with_quality_check<dtype, mode, subsample>(
            aln, reference, arr, params, stats
        );

        if (stats.mod(params.print_every)) {
            stats.aggregate();
            stats.body();
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
            return {size, length, length, 2, 2};
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
    Stats& stats,
    const std::string& name
) : Main(file, fasta, hdf5, mpi, params, stats),
    memspace(__memspace<dtype, _ndims(mode)>(fasta, hdf5, name, mode)) {}

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
    Stats& stats,
    const std::string& name
) {

    switch (params.mode) {

        case Mode::Normal: {

                return std::make_unique<TemplatedMain<float, Mode::Normal>>(
                    file, fasta, hdf5, mpi, params, stats, name
                );

        }

        case Mode::LowMem: {

                return std::make_unique<TemplatedMain<float, Mode::LowMem>>(
                    file, fasta, hdf5, mpi, params, stats, name
                );

        }

        case Mode::Joint: {

                return std::make_unique<TemplatedMain<float, Mode::Joint>>(
                    file, fasta, hdf5, mpi, params, stats, name
                );

        }

        case Mode::Tokenize: {

            throw std::runtime_error("Mode::Tokenize is not a valid input to _get_main.");

        }

    }

}





//
// Stats
//





static inline double _percent(int64_t a, int64_t b) {
    if (b == 0) return 0.0;
    double a_double = static_cast<double>(a);
    double b_double = static_cast<double>(b);
    return a_double / b_double * 100;
}

static inline void _print_header(int64_t files, int64_t aligned, int64_t unaligned, int64_t references, int64_t length) {

    Utils::Line _print_references("References");
    Utils::Line _print_length("Reference length");
    Utils::Line _print_aligned("Aligned reads");
    Utils::Line _print_unaligned("Unaligned reads");

    _print_references.print(references);
    _print_length.print(length);
    _print_aligned.print(aligned);
    _print_unaligned.print(unaligned);

    Utils::divider();
    Utils::cursor_down(4);
    Utils::cursor_down(2);
    Utils::cursor_up(2);

}





//
// Stats
//





Stats::Stats(
    int64_t files,
    int64_t aligned,
    int64_t unaligned,
    int32_t references,
    int32_t length,
    const MPI::Manager& mpi
) 
    : _aligned(aligned), _unaligned(unaligned), _references(references), _length(length), _files(files), _mpi(mpi) {

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


void Stats::file() {
    _curr_file++;
}


void Stats::aggregate() {
    _processed = _mpi.reduce(_processed);
    _skipped   = _mpi.reduce(_skipped);
}


void Stats::header() const {
    if (_mpi.root()) {
        _print_header(_files, _aligned, _unaligned, _references, _length);
    }
}


void Stats::body() const {

    int fields = 4;
    _mpi.up(fields);

    if (_mpi.root()) {
        double processed = _percent(_processed, _aligned + _unaligned);
        double skipped   = _percent(_skipped, _aligned + _unaligned);
        std::string file = std::to_string(_curr_file) + "/" + std::to_string(_files);
        _print_files.print(file);
        _print_processed.print(processed);
        _print_skipped.print(skipped);
        _print_elapsed.print(_mpi.time_str());
    }

    _mpi.divide();
    _mpi.up();

}

bool Stats::mod(int64_t n) const {

    return (_processed % n) == 0;

}





} // namespace cmuts
