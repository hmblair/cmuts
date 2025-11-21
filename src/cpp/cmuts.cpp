#include "cmuts.hpp"


template <typename dtype>
class Stack {
private:

    std::vector<dtype> _data;
    int32_t _ix = -1;

public:

    Stack(int32_t max);
    void push(dtype value);
    dtype pop();

    size_t size() const;
    dtype top() const;
    std::vector<dtype> data() const;
    bool empty() const;

};

template <typename dtype>
Stack<dtype>::Stack(int32_t max)
    : _data(max) {}

template <typename dtype>
void Stack<dtype>::push(dtype value) {
    _ix++;
    if (_ix >= _data.size()) {
        _data.resize(2 * _data.size());
    }
    _data[_ix] = value;
}

template <typename dtype>
dtype Stack<dtype>::pop() {
    if (empty()) {
        throw std::runtime_error("Cannot pop from empty stack.");
    }
    _ix--;
    return _data[_ix + 1];
}

template <typename dtype>
size_t Stack<dtype>::size() const {
    return _ix + 1;
}

template <typename dtype>
std::vector<dtype> Stack<dtype>::data() const {
    return {_data.begin(), _data.begin() + size()};
}

template <typename dtype>
dtype Stack<dtype>::top() const {
    return _data[_ix];
}

template <typename dtype>
bool Stack<dtype>::empty() const {
    return size() == 0;
}





template <typename dtype>
static inline dtype _pow(dtype base, int32_t n) {

    dtype tmp = 1;
    for (int32_t ix = 0; ix < n; ix++) { tmp *= base; }
    return tmp;

}


template <typename dtype>
static inline void _fill_at_indices(
    std::vector<dtype>& vec,
    const std::vector<int32_t>& indices,
    dtype value
) {

    for (const auto& ix: indices) { vec[ix] += value; }

}

#include <iomanip>
template <typename dtype>
void _print_vector(const std::vector<dtype>& vec, int precision = 4) {

    size_t size = vec.size();
    if (size == 0) { return; }

    std::cout << std::fixed << std::setprecision(precision);
    std::cout << "{";
    for (size_t ix = 0; ix < size - 1; ix++) {
        std::cout << vec[ix] << ",";
    }
    std::cout << vec[size - 1] << "}" << std::endl;

}





namespace cmuts {


std::vector<bool> _ignore_str_to_bool(const std::string& ignore) {

    std::vector<bool> bases(BASES, true);
    for (const auto& base: ignore) {
        base_t val = HTS::_from_char(base);
        if (val == IX_UNK) {
            throw std::invalid_argument("Invalid base " + std::string(1, base) + ".");
        }
        bases[val] = false;
    }
    return bases;

}


static inline std::vector<size_t> __dims(const BinaryFASTA& fasta, Mode mode) {

    size_t size   = fasta.size();
    size_t length = fasta.longest();

    switch (mode) {

        case Mode::Normal: {
            return {size, length, N_BASES, N_DELBASES};
        }

        case Mode::Pairwise: {
            return {size, length, length, N_PAIRS, N_PAIRS};
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
DataView<dtype, mode>::DataView(HDF5::Memspace<dtype, _ndims(mode)> memspace) :
    _memspace(memspace) {}


template <typename dtype, Mode mode>
view_t<dtype, _ndims(mode)> DataView<dtype, mode>::view() {

    return _memspace.view(_offset);

}

template <typename dtype, Mode mode>
void DataView<dtype, mode>::update(int32_t offset) {

    _offset = offset;

}

template <typename dtype, Mode mode>
void DataView<dtype, mode>::write(int32_t offset) {

    if (!skip) {
        _memspace.safe_write(offset);
        _memspace.clear();
        skip = true;
    }

}

template <typename dtype, Mode mode>
int32_t DataView<dtype, mode>::size() const {

    return _memspace.size();

}


template <typename dtype>
Data<dtype>::Data(
    BinaryFASTA& fasta,
    HDF5::File& hdf5,
    const std::string& name,
    bool pairwise
) :
    mods(__memspace<dtype, _ndims(Mode::Normal)>(fasta, hdf5, name + "/counts-1d", Mode::Normal))
{

    if (pairwise) {
        pairs = DataView<dtype, Mode::Pairwise>(
            __memspace<dtype, _ndims(Mode::Pairwise)>(fasta, hdf5, name + "/counts-2d", Mode::Pairwise)
        );
    }

}

template <typename dtype>
void Data<dtype>::update(int32_t offset) {

    mods.update(offset);
    if (pairs) { pairs->update(offset); }

}


template <typename dtype>
void Data<dtype>::write(int32_t offset) {

    mods.write(offset);
    if (pairs) { pairs->write(offset); }

}





//
// Helpers for joint
//





template <typename dtype>
static inline void __joint(Data<dtype>& data) {

    auto arr = data.pairs->view();

    #pragma omp parallel for schedule(dynamic)
    for (int32_t ix = data.min; ix < data.tmp.size(); ix++) {

        dtype ix_val   = data.tmp[ix];
        dtype ix_val_c = 1 - ix_val;

        #if defined(__clang__)
            #pragma clang loop vectorize(enable)
        #elif defined(__GNUC__) && !defined(__clang__)
            #pragma GCC ivdep
        #endif
        for (int32_t jx = data.min; jx < data.tmp.size(); jx++) {

            dtype jx_val   = data.tmp[jx];
            dtype jx_val_c = 1 - jx_val;

            arr(ix, jx, NOMOD, NOMOD) += ix_val_c * jx_val_c;
            arr(ix, jx, NOMOD, MOD)   += ix_val_c * jx_val;
            arr(ix, jx, MOD, NOMOD)   += ix_val * jx_val_c;
            arr(ix, jx, MOD, MOD)     += ix_val * jx_val;

        }

    }

}





//
// For ambiguous deletions
//





template <typename dtype>
dtype _deletion_prob(
    const std::vector<dtype>& mutations,
    const std::vector<int32_t>& indices,
    dtype penalty
) {

    auto n = static_cast<int32_t>(indices.size() - 1);
    dtype tmp = _pow(penalty, n);

    for (const auto& ix: indices) { tmp *= mutations[ix]; }
    return tmp;

}


struct DeletionData {

    int32_t M = 0;
    std::vector<int32_t> ix;

};


int32_t _get_ambiguous_end(
    int32_t start,
    int32_t end,
    const seq_t& sequence
) {

    auto size = static_cast<int32_t>(sequence.size());
    if (end >= size) { return size; }

    int32_t M = start;
    int32_t N = end;

    while (M < N && N < size) {
        if (sequence[M] == sequence[N]) { N++; }
        M++;
    }

    // This is the position of the base _after_ the final base that can be deleted.
    return N;

}


template <typename dtype>
static inline dtype _total_muts(
    Data<dtype>& data,
    int32_t ix
) {

    dtype total = 0;
    dtype diag  = 0;

    auto _view = data.mods.view();

    // From A
    diag  += _view(ix, IX_A, IX_A);
    total += _view(ix, IX_A, IX_C);
    total += _view(ix, IX_A, IX_G);
    total += _view(ix, IX_A, IX_T);
    // From C
    total += _view(ix, IX_C, IX_A);
    diag  += _view(ix, IX_C, IX_C);
    total += _view(ix, IX_C, IX_G);
    total += _view(ix, IX_C, IX_T);
    // From G
    total += _view(ix, IX_G, IX_A);
    total += _view(ix, IX_G, IX_C);
    diag  += _view(ix, IX_G, IX_G);
    total += _view(ix, IX_G, IX_T);
    // From T/U
    total += _view(ix, IX_T, IX_A);
    total += _view(ix, IX_T, IX_C);
    total += _view(ix, IX_T, IX_G);
    diag  += _view(ix, IX_T, IX_T);

    dtype coverage = total + diag;

    if (coverage > 0) {
        return total / coverage;
    } else {
        return total;
    }

}


template <typename dtype>
static inline void _normalize(std::vector<dtype>& arr, dtype sum = 1) {

    dtype total = 0;
    for (const dtype& val : arr) { total += val; }

    if (total == 0) {
        for (dtype& val : arr) { val = sum / arr.size(); }
    } else {
        total /= sum;
        for (dtype& val : arr) { val /= total; }
    }

}


template <typename dtype>
static inline std::vector<dtype> _spread_weights(
    int32_t start,
    int32_t end,
    Data<dtype>& data,
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

            for (int32_t ix = 0; ix < length; ix++) {
                weights[ix] = mask / length;
            }
            break;

        }

        case Spread::MutationInformed: {

            for (int32_t ix = start; ix < end; ix++) {
                weights[ix - start] = _total_muts<dtype>(data, ix) * mask;
            }

            break;

        }

    }

    return weights;

}


template <typename dtype>
std::vector<dtype> _enumerate_deletions(
    const seq_t& reference,
    const std::vector<dtype>& mutations,
    dtype penalty,
    dtype mask,
    int32_t start,
    int32_t end
) {

    // Start of the search region

    int32_t M = start;

    // End of the search region (and the current deletion)

    int32_t N = end;

    // Core stack

    int32_t max = 32;
    Stack<DeletionData> stack(max);

    DeletionData init;
    init.M = M;
    init.ix = {N - 1 - start};
    stack.push(init);

    // Final deletion profile

    std::vector<dtype> out(mutations.size(), 0);

    // Is the current value of M on the boundary of the current deletion

    bool contig = true;

    // Did we push to the stack this loop

    bool pushed = false;

    // Keep track of the average number of mutations
    // for normalization purposes

    dtype total = 0;
    dtype nmuts = 0;

    while (!stack.empty()) {

        pushed = false;

        while (M < N && N < reference.size()) {

            if (reference[M] == reference[N]) {

                M++;
                N++;

                DeletionData data;
                data.M = M;
                data.ix = stack.top().ix;

                // The top of the index stack gets replaced

                if (contig) {

                    data.ix.back() = N - 1 - start;

                }

                // The top of the index stack gets split in two

                else {

                    data.ix.back() = stack.top().M - start;
                    data.ix.push_back(N - 1 - start);

                }

                stack.push(data);
                pushed = true;
                contig = true;

                break;

            }

            // Advance M

            M++;
            contig = false;

        }

        if (!pushed) {

            // Pop and accumulate

            auto data = stack.pop();
            M = data.M - 1;
            N--;

            float prob = _deletion_prob<dtype>(mutations, data.ix, penalty);
            _fill_at_indices<dtype>(out, data.ix, prob);

            nmuts += data.ix.size() * prob;
            total += prob;

            // Advance M

            M++;
            contig = false;

        }

    }

    if (total > 0) { nmuts /= total; }
    else { nmuts = 1; }

    _normalize<dtype>(out, mask * nmuts);
    return out;

}





//
// Mutation counting
//





template <typename dtype>
static inline void __match_core(
    Data<dtype>& data,
    int32_t rpos,
    base_t rbase,
    dtype mask,
    const Params& params
) {

    data.mods.view()(rpos, rbase, rbase) += mask;

};


template <typename dtype>
static inline void __mismatch_core(
    Data<dtype>& data,
    int32_t rpos,
    base_t qbase,
    base_t rbase,
    dtype mask,
    const Params& params
) {

    if (!params.bases[rbase]) { mask = 0; }

    data.mods.view()(rpos, rbase, qbase) += mask;
    if (params.pairwise) { data.tmp[rpos] += mask; }

};

template <typename dtype>
static inline void __ins_core(
    Data<dtype>& data,
    int32_t rpos,
    base_t qbase,
    dtype mask,
    const Params& params
) {

    data.mods.view()(rpos, qbase, IX_INS) += mask;
    if (params.pairwise) { data.tmp[rpos] += mask; }

};

template <typename dtype>
static inline void __del_core(
    Data<dtype>& data,
    int32_t rpos,
    base_t rbase,
    dtype mask,
    const Params& params
) {

    if (!params.bases[rbase]) { mask = 0; }

    data.mods.view()(rpos, rbase, IX_DEL) += mask;
    if (params.pairwise) { data.tmp[rpos] += mask; }

};


template <typename dtype>
static inline void __term(
    Data<dtype>& data,
    int32_t rpos,
    base_t rbase
) {

    data.mods.view()(rpos, rbase, IX_TERM) += 1;

};


template <
    typename dtype,
    bool CONSUMES_RPOS,
    bool CONSUMES_QPOS
>
static inline void __match(
    Data<dtype>& data,
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

        dtype value = mask[qpos];
        if (params.no_filter_matches) { value = 1; }

        __match_core<dtype>(data, rpos, rbase, value, params);

    }

}


template <typename dtype>
static inline void __mismatch(
    Data<dtype>& data,
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
    dtype value  = mask[qpos];

    if (params.mismatches && last - rpos >= params.collapse && qbase != IX_UNK) {
        __mismatch_core<dtype>(data, rpos, qbase, rbase, value, params);
        last = rpos;
    } else {
        if (params.no_filter_matches) { value = 1; }
        __match_core<dtype>(data, rpos, rbase, value, params);
    }

    __match<dtype, true, true>(data, op, rpos, qpos, reference, mask, params);

}


template <typename dtype>
static inline void __ins(
    Data<dtype>& data,
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
        !params.insertions                    ||
        op.length() > params.max_indel_length ||
        rpos > reference.size()               ||
        rpos <= 0
    ) {
        qpos -= op.length();
        return;
    }

    qpos--;
    base_t qbase = op.last();

    dtype value = mask[qpos];
    if (params.no_filter_insertions) { value = 1; }

    if (last - rpos >= params.collapse && qbase != IX_UNK) {
        __ins_core<dtype>(data, rpos - 1, qbase, value, params);
        last = rpos;
    }

    qpos -= (op.length() - 1);

}

template <typename dtype>
static inline void __del_no_ambig(
    Data<dtype>& data,
    HTS::CIGAR_op& op,
    int32_t& rpos,
    int32_t& qpos,
    int32_t& last,
    const seq_t& reference,
    const std::vector<dtype>& mask,
    const Params& params
) {

    op.advance();
    rpos--;
    last = rpos;

    base_t rbase = reference[rpos];

    dtype value = mask[qpos];
    if (params.no_filter_deletions) { value = 1; }

    __del_core<dtype>(data, rpos, rbase, value, params);

    __match<dtype, true, false>(data, op, rpos, qpos, reference, mask, params);

}


template <typename dtype>
static inline void __del_w_ambig(
    Data<dtype>& data,
    HTS::CIGAR_op& op,
    int32_t& rpos,
    int32_t& qpos,
    int32_t& last,
    const seq_t& reference,
    const std::vector<dtype>& mask,
    const Params& params
) {

    // Start and end indices of the deletion

    int32_t start = rpos - op.length();
    int32_t end   = _get_ambiguous_end(start, rpos, reference);

    // RT repeat probability

    dtype penalty = 1 / 0.5;

    dtype value = mask[qpos];
    if (params.no_filter_deletions) { value = 1; }

    std::vector<dtype> mutations = _spread_weights(
        start, end, data, value, params.spread
    );

    std::vector<dtype> deletions = _enumerate_deletions<dtype>(
        reference, mutations, penalty, value, start, rpos
    );

    for (int32_t ix = end - 1; ix >= start; ix--) {
        base_t rbase = reference[ix];
        __del_core<dtype>(data, ix, rbase, deletions[ix - start], params);
    }

    op.advance();
    rpos--;
    last = rpos;

    __match<dtype, true, false>(data, op, rpos, qpos, reference, mask, params);

}


template <typename dtype>
static inline void __del(
    Data<dtype>& data,
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
        __match<dtype, true, false>(data, op, rpos, qpos, reference, mask, params);
        return;
    }

    if (params.ambiguous) {
        __del_w_ambig<dtype>(data, op, rpos, qpos, last, reference, mask, params); 
    } else {
        __del_no_ambig<dtype>(data, op, rpos, qpos, last, reference, mask, params);
    }

}


template <typename dtype>
static inline void __count(
    HTS::Alignment& aln,
    const seq_t& reference,
    Data<dtype>& data,
    const Params& params,
    Stats& stats
) {

    HTS::CIGAR& cigar       = aln.cigar;
    const HTS::PHRED& phred = aln.phred;

    // Initial reference and query positions

    int32_t rpos = aln.offset + cigar.rlength();
    int32_t qpos = aln.length;

    if (rpos > reference.size()) {
        __throw_and_log(
            _LOG_FILE,
            "Reference position " + std::to_string(rpos) + 
            " exceeds reference sequence length " + std::to_string(reference.size())
        );
    }

    // Fill the tmp modification array if necessary

    if (params.pairwise) { data.tmp.assign(rpos, 0); }

    // Position of the most recent (3'-wise) mutation

    int32_t last = rpos + params.collapse;

    // PHRED quality mask

    HTS::BaseMask<dtype> mask = phred.mask<dtype>(
        params.min_phred,
        params.quality_window
    );

    // Update number of low-quality bases

    stats.update_bases(mask.mask.size() - mask.good, mask.mask.size());

    // Loop backwards (3' -> 5') along the CIGAR

    for (auto& op : std::ranges::reverse_view(cigar)) {

        switch (op.type()) {

            case HTS::CIGAR_t::MATCH: {

                __match<dtype, true, true>(data, op, rpos, qpos, reference, mask.mask, params);
                break;

            }

            case HTS::CIGAR_t::MISMATCH: {

                __mismatch<dtype>(data, op, rpos, qpos, last, reference, mask.mask, params);
                break;

            }

            case HTS::CIGAR_t::DEL: {

                __del<dtype>(data, op, rpos, qpos, last, reference, mask.mask, params);
                break;

            }

            case HTS::CIGAR_t::INS: {

                __ins<dtype>(data, op, rpos, qpos, last, reference, mask.mask, params);
                break;

            }

            case HTS::CIGAR_t::TERM: {

                __term<dtype>(data, rpos, reference[rpos]);
                break;

            }

            case HTS::CIGAR_t::SOFT:
            case HTS::CIGAR_t::SKIP: {

                qpos -= op.length();
                break;

            }

            case HTS::CIGAR_t::HARD:
            case HTS::CIGAR_t::PAD:
            case HTS::CIGAR_t::BACK: {

                break;

            }

            case HTS::CIGAR_t::UNKNOWN: {

                __throw_and_log(
                    _LOG_FILE,
                    "Unknown CIGAR operation encountered at reference position "
                    + std::to_string(rpos) + "."
                );

            }

        }

    }

    // Store where the read stopped

    if (params.pairwise) { data.min = rpos; }

}


static inline bool __check_sense(
    const HTS::Alignment& aln,
    const Params& params
) {

    return (params.forward || aln.reversed) && (params.reverse || !aln.reversed);

}


static inline bool __check_mapq(
    const HTS::Alignment& aln,
    const Params& params
) {

    return (aln.mapq >= params.min_mapq) && (aln.mapq < MISSING_MAPQ || params.min_mapq == 0);

}


static inline bool __check_primary(
    const HTS::Alignment& aln,
    const Params& params
) {

    return (aln.primary || params.secondary);

}


static inline bool __check_length(
    const HTS::Alignment& aln,
    const Params& params
) {

    return (aln.length >= params.min_length) && (aln.length <= params.max_length);

}


static inline bool __check_hamming(
    const HTS::Alignment& aln,
    const Params& params
) {

    return aln.cigar.hamming() <= params.max_hamming;

}


static inline bool __check_all(
    const HTS::Alignment& aln,
    const Params& params
) {

    return (
        aln.aligned                  &&
        __check_primary(aln, params) &&
        __check_sense(aln, params)   &&
        __check_mapq(aln, params)    &&
        __check_length(aln, params)  &&
        __check_hamming(aln, params)
    );

}


template <typename dtype>
static inline void __count_with_quality_check(
    HTS::Alignment& aln,
    const seq_t& reference,
    Data<dtype>& data,
    const Params& params,
    Stats& stats
) {

    if (!__check_all(aln, params)) {
        stats.skipped();
    } else {
        __count<dtype>(aln, reference, data, params, stats);
        if (params.pairwise) { __joint<dtype>(data); }
        stats.processed();
    }

}


template <typename dtype>
static inline void __count_reference(
    HTS::Iterator& iter,
    const seq_t& reference,
    Data<dtype>& data,
    const Params& params,
    Stats& stats
) {

    int32_t count = 0;
    while (!iter.end() && count < params.downsample) {

        HTS::Alignment aln = iter.next();
        __count_with_quality_check<dtype>(aln, reference, data, params, stats);

        stats.print();
        count++;

        data.mods.skip = false;
        if (params.pairwise) { data.pairs->skip = false; }

    }

}





//
// Main
//





template <typename dtype>
static inline void __process(
    HTS::File& file,
    BinaryFASTA& fasta,
    HDF5::File& hdf5,
    const Params& params,
    Stats& stats,
    Data<dtype>& data,
    int32_t min
) {

    int32_t max = std::min(
        min + data.mods.size(),
        file.size()
    );

    for (int32_t ix = min; ix < max; ix++) {

        // Reads are processed sequentially if built in serial mode, so
        // no need to seek at beginning of new chunk

        #ifdef MPI_BUILD
        bool seek = (ix % hdf5.chunk_size() == 0);
        #else
        bool seek = (ix == 0);
        #endif

        std::shared_ptr<HTS::Iterator> iter = file.get(ix, seek);
        if (iter->end()) { continue; }

        seq_t reference = fasta.sequence(ix);
        data.update(ix - min);
        __count_reference<dtype>(*iter, reference, data, params, stats);

    }

    data.write(min);

}


template <typename dtype>
void run(
    HTS::File& file,
    BinaryFASTA& fasta,
    HDF5::File& hdf5,
    const MPI::Manager& mpi,
    const Params& params,
    Stats& stats,
    const std::string& name
) {

    stats.body();
    MPI::Chunk chunk = mpi.chunk(hdf5.chunk_size(), file.size());

    Data<dtype> data(fasta, hdf5, name, params.pairwise);

    for (int32_t min = chunk.low; min < chunk.high; min += chunk.step) {

        __process<dtype>(
            file,
            fasta,
            hdf5,
            params,
            stats,
            data,
            min
        );

    }

}


template void run<float>(
    HTS::File& file,
    BinaryFASTA& fasta,
    HDF5::File& hdf5,
    const MPI::Manager& mpi,
    const Params& params,
    Stats& stats,
    const std::string& name
);


Spread spread(bool uniform, bool none) {

    if (uniform && none) {
        __throw_and_log(
            _LOG_FILE,
            "--uniform-spread and --no-spread are mutually exclusive."
        );
    }

    if (uniform) { return Spread::Uniform;          }
    if (none)    { return Spread::None;             }
    else         { return Spread::MutationInformed; }

}





//
// Stats
//





static inline double _percent(int64_t a, int64_t b) {

    double max = 100;
    if (b == 0) return max;
    return static_cast<double>(a) / static_cast<double>(b) * max;

}

static inline void _print_header(
    int64_t files,
    int64_t aligned,
    int64_t unaligned,
    int64_t references,
    int64_t length
) {

    Utils::Line _print_references("References");
    Utils::Line _print_length("Reference length");
    Utils::Line _print_aligned("Aligned reads");
    Utils::Line _print_unaligned("Unaligned reads");

    _print_references.print(references);
    _print_length.print(length);
    _print_aligned.print(aligned);
    _print_unaligned.print(unaligned);

    Utils::divider();
    Utils::cursor_down(5);
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
    double print_every,
    const MPI::Manager& mpi
) : 
    _aligned(aligned),
    _unaligned(unaligned),
    _references(references),
    _length(length),
    _files(files),
    _print_every(print_every),
    _mpi(mpi) {}


void Stats::update_bases(int64_t skipped, int64_t total) {

    _bases_skipped   += skipped;
    _bases_processed += total;

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

    _bases_processed = _mpi.reduce(_bases_processed);
    _bases_skipped   = _mpi.reduce(_bases_skipped);

}


void Stats::header() const {
    if (_mpi.root()) {
        _print_header(_files, _aligned, _unaligned, _references, _length);
    }
}


void Stats::body() {

    int fields = 5;
    _mpi.up(fields);

    if (_mpi.root()) {

        double processed     = _percent(_processed, _aligned);
        double skipped       = _percent(_skipped, _processed);
        double bases_skipped = _percent(_bases_skipped, _bases_processed);

        std::string file = std::to_string(_curr_file) + "/" + std::to_string(_files);

        _print_files.print(file);
        _print_processed.print(processed);
        _print_skipped.print(skipped);
        _print_bases_skipped.print(bases_skipped);
        _print_elapsed.print(_mpi.time_str());

    }

    _mpi.divide();
    _mpi.up();

    _last_print = _mpi.time();

}

double Stats::elapsed() const {

    return _mpi.time() - _last_print;

}

void Stats::print() {

    if (elapsed() > _print_every) {
        aggregate();
        body();
    }

}





} // namespace cmuts
