#include "tests.hpp"



//
// For generation
//



static inline std::mt19937 _init_gen() {
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    return std::mt19937(seed);
}

static inline int _random_int(std::mt19937& gen, int low, int high) {

    std::uniform_int_distribution<> dist(low, high);
    return dist(gen);

}

static inline int _random_int_excluding(std::mt19937& gen, int low, int high, int n) {

    int rand_val = _random_int(gen, low, high - 1);
    if (rand_val >= n) {
        ++rand_val;
    }
    return rand_val;

}

template <typename dtype>
static inline dtype _sample_from_vector(std::mt19937& gen, const std::vector<dtype>& vals) {

    int ix = _random_int(gen, 0, vals.size() - 1);
    return vals[ix];

}

static inline base_t _random_base(std::mt19937& gen) {

    return _random_int(gen, 0, 3);

}

static inline base_t _random_base_excluding(std::mt19937& gen, base_t n) {

    return _random_int_excluding(gen, 0, 3, n);

}

static inline seq_t _random_sequence(
    hts_pos_t length,
    std::mt19937& gen
) {

    seq_t sequence(length);
    for (hts_pos_t ix = 0; ix < length; ix++) {
        sequence[ix] = _random_base(gen);
    }

    return sequence;

}

static inline int _random_length(std::mt19937& gen, size_t max) {

    return _random_int(gen, 1, max);

}

static inline int _random_mapping_quality(std::mt19937& gen) {

    return _random_int(gen, 0, MAX_MAPQ);

}

static inline std::vector<int> _random_mapping_qualities(std::mt19937& gen, size_t count) {

    std::vector<int> _mapping_qualities(count);
    for (auto& val : _mapping_qualities) {
        val = _random_mapping_quality(gen);
    }

    return _mapping_qualities;

}

static inline uint8_t _random_phred_quality(std::mt19937& gen) {

    return _random_int(gen, 0, MAX_PHRED);

}

static inline std::vector<uint8_t> _random_phred_qualities(std::mt19937& gen, size_t count) {

    std::vector<uint8_t> _phred_qualities(count);
    for (auto& val : _phred_qualities) {
        val = _random_phred_quality(gen);
    }

    return _phred_qualities;

}

static inline HTS::CIGAR_op _random_cigar_op(
    std::mt19937& gen,
    hts_pos_t match_remaining,
    hts_pos_t ins_remaining,
    hts_pos_t del_remaining,
    HTS::CIGAR_t prev
) {

    HTS::CIGAR_t type;
    hts_pos_t length;

    if (match_remaining == 0 && del_remaining == 0) {
        type = HTS::CIGAR_t::INS;
        length = ins_remaining;
        return HTS::CIGAR_op(type, length);
    }

    if (match_remaining == 0 && ins_remaining == 0) {
        type = HTS::CIGAR_t::DEL;
        length = del_remaining;
        return HTS::CIGAR_op(type, length);
    }

    std::vector<HTS::CIGAR_t> samples;
    if (match_remaining > 0) {
        if (prev != HTS::CIGAR_t::MATCH) {
            samples.push_back(HTS::CIGAR_t::MATCH);
        }
        if (prev != HTS::CIGAR_t::MISMATCH) {
            samples.push_back(HTS::CIGAR_t::MISMATCH);
        }
    }
    if (ins_remaining > 0) {
        if (prev != HTS::CIGAR_t::INS) {
            samples.push_back(HTS::CIGAR_t::INS);
        }
    }
    if (del_remaining > 0) {
        if (prev != HTS::CIGAR_t::DEL) {
            samples.push_back(HTS::CIGAR_t::DEL);
        }
    }

    type = _sample_from_vector<HTS::CIGAR_t>(gen, samples);
    if (type == HTS::CIGAR_t::MATCH || type == HTS::CIGAR_t::MISMATCH) {
        length = _random_length(gen, match_remaining);
    } else if (type == HTS::CIGAR_t::INS) {
        length = _random_length(gen, ins_remaining);
    } else if (type == HTS::CIGAR_t::DEL) {
        length = _random_length(gen, del_remaining);
    } else {
        throw std::runtime_error("Error in CIGAR type.");
    }

    return HTS::CIGAR_op(type, length);

}




template <typename dtype>
void _add_matches(
    const seq_t& reference,
    seq_t& query,
    hts_pos_t& qpos,
    hts_pos_t& rpos,
    view_t<dtype, 4> arr,
    int quality,
    int min_quality,
    const HTS::PHRED& phred,
    size_t window,
    HTS::CIGAR_op& op,
    std::mt19937& gen
) {


    while (op.advance()) {

        bool _mapping_mask = (quality >= min_quality);
        bool _phred_mask = phred.check(qpos, min_quality, window);
        dtype mask = static_cast<dtype>(_mapping_mask && _phred_mask);

        base_t rbase = reference[rpos];
        query.push_back(rbase);
        arr(rpos, rbase, rbase) += mask;

        qpos++;
        rpos++;

    }

}

template <typename dtype>
void _add_mismatches(
    const seq_t& reference,
    seq_t& query,
    hts_pos_t& qpos,
    hts_pos_t& rpos,
    view_t<dtype, 4> arr,
    int quality,
    int min_quality,
    const HTS::PHRED& phred,
    size_t window,
    HTS::CIGAR_op& op,
    std::mt19937& gen
) {

    while (op.advance(-1)) {

        bool _mapping_mask = (quality >= min_quality);
        bool _phred_mask = phred.check(qpos, min_quality, window);
        dtype mask = static_cast<dtype>(_mapping_mask && _phred_mask);

        base_t rbase = reference[rpos];
        base_t qbase = _random_base_excluding(gen, rbase);
        query.push_back(qbase);

        arr(rpos, rbase, rbase) += mask;

        qpos++;
        rpos++;
    }

    base_t rbase = reference[rpos];
    base_t qbase = _random_base_excluding(gen, rbase);
    query.push_back(qbase);

    bool _mapping_mask = (quality >= min_quality);
    bool _phred_mask = phred.check(qpos, min_quality, window);
    dtype mask = static_cast<dtype>(_mapping_mask && _phred_mask);

    arr(rpos, rbase, qbase) += mask;

    qpos++;
    rpos++;

}

template <typename dtype>
void _add_deletions(
    const seq_t& reference,
    seq_t& query,
    hts_pos_t& qpos,
    hts_pos_t& rpos,
    view_t<dtype, 4> arr,
    int quality,
    int min_quality,
    const HTS::PHRED& phred,
    size_t window,
    HTS::CIGAR_op& op,
    std::mt19937& gen,
    bool _matches_instead
) {

    bool _mapping_mask = (quality >= min_quality);
    bool _phred_mask = phred.check(qpos, min_quality, window);
    dtype mask = static_cast<dtype>(_mapping_mask && _phred_mask);

    while (op.advance(-1)) {

        base_t rbase = reference[rpos];
        arr(rpos, rbase, rbase) += mask;
        rpos++;

    }

    base_t rbase = reference[rpos];
    if (_matches_instead) {
        arr(rpos, rbase, rbase) += mask;
    } else {
        arr(rpos, rbase, IX_DEL) += mask;
    }

    rpos++;

}

template <typename dtype>
void _add_insertions(
    const seq_t& reference,
    seq_t& query,
    hts_pos_t& qpos,
    hts_pos_t& rpos,
    view_t<dtype, 4> arr,
    int quality,
    int min_quality,
    const HTS::PHRED& phred,
    size_t window,
    HTS::CIGAR_op& op,
    std::mt19937& gen
) {

    while (op.advance(-1)) {

        base_t qbase = _random_base(gen);
        query.push_back(qbase);
        qpos++;

    }

    bool _mapping_mask = (quality >= min_quality);
    bool _phred_mask = phred.check(qpos, min_quality, window);
    dtype mask = static_cast<dtype>(_mapping_mask && _phred_mask);

    base_t qbase = _random_base(gen);
    query.push_back(qbase);
    if (rpos < reference.size()) {
        arr(rpos, qbase, IX_INS) += mask;
    }

    qpos++;

}

template <typename dtype>
struct _MUTANT {
    seq_t query;
    HTS::CIGAR cigar;
    std::string phred;
};

template <typename dtype>
_MUTANT<dtype> _get_mutant(
    const seq_t& reference,
    view_t<dtype, 4> arr,
    int quality,
    int min_quality,
    hts_pos_t min_length,
    int max_indel_length,
    size_t window,
    std::mt19937& gen
) {

    seq_t query;
    hts_pos_t qpos = 0, rpos = 0;
    HTS::CIGAR_op curr, prev;
    HTS::CIGAR cigar;

    hts_pos_t ins_length = _random_int(gen, 0, reference.size());
    hts_pos_t del_length = _random_int(gen, 0, reference.size() - 1);
    hts_pos_t match_length = reference.size() - del_length;
    hts_pos_t query_length = match_length + ins_length;

    if (query_length < min_length) { quality = -1; }

    hts_pos_t curr_ins = 0, curr_del = 0, curr_match = 0;

    std::vector<uint8_t> __phred = _random_phred_qualities(gen, query_length);
    HTS::PHRED phred(__phred);

    while (curr_match < match_length || curr_ins < ins_length || curr_del < del_length) {

        hts_pos_t match_remaining = match_length - curr_match;
        hts_pos_t ins_remaining = ins_length - curr_ins;
        hts_pos_t del_remaining = del_length - curr_del;
        curr = _random_cigar_op(gen, match_remaining, ins_remaining, del_remaining, prev.type());

        switch (curr.type()) {

            case HTS::CIGAR_t::MATCH: {

                _add_matches<dtype>(
                    reference,
                    query,
                    qpos,
                    rpos,
                    arr,
                    quality,
                    min_quality,
                    phred,
                    window,
                    curr,
                    gen
                );
                prev = curr;
                cigar.extend(curr);
                curr_match += curr.length();
                break;

            }

            case HTS::CIGAR_t::MISMATCH: {

                _add_mismatches<dtype>(
                    reference,
                    query,
                    qpos,
                    rpos,
                    arr,
                    quality,
                    min_quality,
                    phred,
                    window,
                    curr,
                    gen
                );
                prev = curr;
                // Convert the mismatch to a match as they are the same in the CIGAR str
                cigar.extend(curr.match());
                curr_match += curr.length();
                break;

            }

            case HTS::CIGAR_t::DEL: {

                bool _matches_instead = false;
                if (curr.length() > max_indel_length) { _matches_instead = true; }

                _add_deletions<dtype>(
                    reference,
                    query,
                    qpos,
                    rpos,
                    arr,
                    quality,
                    min_quality,
                    phred,
                    window,
                    curr,
                    gen,
                    _matches_instead
                );
                prev = curr;
                cigar.extend(curr);
                curr_del += curr.length();
                break;

            }

            case HTS::CIGAR_t::INS: {

                int __quality = quality;
                if (curr.length() > max_indel_length) { __quality = -1; }
                _add_insertions<dtype>(
                    reference,
                    query,
                    qpos,
                    rpos,
                    arr,
                    __quality,
                    min_quality,
                    phred,
                    window,
                    curr,
                    gen
                );
                prev = curr;
                cigar.extend(curr);
                curr_ins += curr.length();
                break;

            }

            default: {
                throw std::runtime_error("Error generating test case.");
            }

        }

    }

    _MUTANT<dtype> mutant;
    mutant.query = query;
    mutant.cigar = cigar;
    mutant.phred = phred.str();

    return mutant;

}


void _write_to_sam(
    std::vector<std::string> references,
    std::vector<std::vector<std::string>> queries,
    std::vector<std::vector<int>> mapping_qualities,
    std::vector<std::vector<std::string>> phreds,
    std::vector<std::vector<std::string>> cigars,
    std::string filename
) {

    std::ofstream f(filename);

    // SAM header
    f << "@HD\tVN:1.6\tSO:coordinate\n";
    for (size_t n = 0; n < references.size(); ++n) {
        f << "@SQ\tSN:ref" << n << "\tLN:" << references[n].length() << "\n";
    }

    // SAM body
    for (size_t n = 0; n < references.size(); ++n) {
        for (size_t i = 0; i < queries[n].size(); ++i) {

            std::string read_name = "read" + std::to_string(i + 1);
            int flag = 0;
            std::string rname = "ref" + std::to_string(n);
            int pos = 1;
            int mapq = mapping_qualities[n][i];
            std::string cigar = cigars[n][i];
            std::string rnext = "*";
            int pnext = 0;
            int tlen = 0;
            std::string qual = phreds[n][i];
            std::string query = queries[n][i];
            f << read_name << "\t" << flag << "\t" << rname << "\t" << pos << "\t" << mapq << "\t" << cigar << "\t"
              << rnext << "\t" << pnext << "\t" << tlen << "\t" << query << "\t" << qual << "\n";
        }
    }
}

static inline std::string _path(const std::string& name) {
    std::filesystem::path path(name);
    std::string stem = path.stem().string();
    std::string parent = path.parent_path().string();
    if (parent.empty()) {
        return stem;
    }
    return parent + "/" + stem;
}

void __run_tests(
    const MPI::Manager& mpi,
    size_t nrefs,
    size_t nqueries,
    hts_pos_t length,
    hts_pos_t min_length,
    int min_quality,
    int max_indel_length,
    size_t window,
    std::string out_sam,
    std::string out_fasta,
    std::string out_h5
) {

    std::mt19937 gen = _init_gen();

    std::optional<HDF5::File> _hdf5;
    try {
        _hdf5.emplace(out_h5, HDF5::RWC, mpi, nrefs, 3);
    } catch (const std::exception& e) {
        mpi.err() << "Error: " << e.what() << "\n";
    }
    HDF5::File hdf5 = std::move(_hdf5.value());


    std::vector<size_t> dims = {nrefs, static_cast<size_t>(length), 4, 6};
    HDF5::Memspace memspace = hdf5.memspace<float, 4>(dims, _path(out_sam));


    std::vector<std::string> references;
    std::vector<std::vector<std::string>> queries;
    std::vector<std::vector<int>> mapping_qualities;
    std::vector<std::vector<std::string>> phreds;
    std::vector<std::vector<std::string>> cigars;

    for (size_t ix = 0; ix < nrefs; ix++) {

        view_t<float, 4> arr = memspace.view(ix);

        seq_t _reference = _random_sequence(length, gen);
        references.push_back(HTS::_to_str(_reference));

        std::vector<int> _mapping_qualties = _random_mapping_qualities(gen, nqueries);
        mapping_qualities.push_back(_mapping_qualties);

        std::vector<std::string> _queries;
        std::vector<std::string> _cigars;
        std::vector<std::string> _phreds;
        for (size_t jx = 0; jx < nqueries; jx++) {

            _MUTANT<float> mut = _get_mutant<float>(_reference, arr, _mapping_qualties[jx], min_quality, min_length, max_indel_length, window, gen);

            _queries.push_back(HTS::_to_str(mut.query));
            _cigars.push_back(mut.cigar.str());
            _phreds.push_back(mut.phred);

        }

        queries.push_back(_queries);
        cigars.push_back(_cigars);
        phreds.push_back(_phreds);

    }

    HTS::_default_write_to_fasta(out_fasta, references);
    _write_to_sam(references, queries, mapping_qualities, phreds, cigars, out_sam);
    memspace.safe_write(0);

}

int main(int argc, char** argv) {

    // Initialize MPI processes
    MPI::Manager mpi(argc, argv);
    // Initialize HDF5 manager
    HDF5::Manager _hdf5_manager;
    // Disable native HTS logging as it does not work well in the
    // multi-threaded environment
    HTS::__disable_logging();

    // Parse command line arguments
    testsProgram opt;
    try {
        opt.parse(argc, argv);
    } catch (const std::exception& err) {
        mpi.err() << "Error: " << err.what() << "\n";
        return EXIT_FAILURE;
    }

    try {
        __run_tests(mpi, opt.references, opt.queries, opt.length, opt.min_length, opt.min_quality, opt.max_indel_length, opt.window, opt.out_sam, opt.out_fasta, opt.out_h5);
        } catch (const std::exception& err) {
        mpi.err() << "Error: " << err.what() << "\n";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;

}
