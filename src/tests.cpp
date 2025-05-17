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
    if (rand_val >= n) { ++rand_val; }
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

static inline int _random_mapq(std::mt19937& gen) {

    return _random_int(gen, 0, MAX_MAPQ);

}

static inline uint8_t _random_phred_quality(std::mt19937& gen) {

    return _random_int(gen, 0, MAX_PHRED);

}

static inline TinyHTS::PHRED _random_phred(std::mt19937& gen, size_t count) {

    std::vector<uint8_t> _phred(count);
    for (auto& val : _phred) {
        val = _random_phred_quality(gen);
    }

    TinyHTS::PHRED phred(_phred);
    return phred;

}

static inline TinyHTS::CIGAR_op _random_cigar_op(
    std::mt19937& gen,
    hts_pos_t match,
    hts_pos_t ins,
    hts_pos_t del,
    TinyHTS::CIGAR_t prev
) {

    TinyHTS::CIGAR_t type;
    hts_pos_t length;

    if (match == 0) {
        if (ins == 0) {
            return TinyHTS::CIGAR_op(TinyHTS::CIGAR_t::DEL, del);
        }
        if (del == 0) {
            return TinyHTS::CIGAR_op(TinyHTS::CIGAR_t::INS, ins);
        }
    }

    std::vector<TinyHTS::CIGAR_t> samples;
    if (match > 0) {
        if (prev != TinyHTS::CIGAR_t::MATCH) {
            samples.push_back(TinyHTS::CIGAR_t::MATCH);
        }
        if (prev != TinyHTS::CIGAR_t::MISMATCH) {
            samples.push_back(TinyHTS::CIGAR_t::MISMATCH);
        }
    }
    if (ins > 0) {
        if (prev != TinyHTS::CIGAR_t::INS) {
            samples.push_back(TinyHTS::CIGAR_t::INS);
        }
    }
    if (del > 0) {
        if (prev != TinyHTS::CIGAR_t::DEL) {
            samples.push_back(TinyHTS::CIGAR_t::DEL);
        }
    }

    type = _sample_from_vector<TinyHTS::CIGAR_t>(gen, samples);

    switch (type) {

        case TinyHTS::CIGAR_t::MATCH:
        case TinyHTS::CIGAR_t::MISMATCH: {

            length = _random_length(gen, match);
            break;

        }

        case TinyHTS::CIGAR_t::INS:      {

            length = _random_length(gen, ins);
            break;

        }

        case TinyHTS::CIGAR_t::DEL:      {

            length = _random_length(gen, del);
            break;

        }

        default: {

            throw std::runtime_error("Error in CIGAR type.");

        }

    }

    return TinyHTS::CIGAR_op(type, length);

}

template <typename dtype>
class Record {
private:

    hts_pos_t _match;
    hts_pos_t _ins;
    hts_pos_t _del;

    hts_pos_t _rpos;
    hts_pos_t _qpos;
    hts_pos_t _last_mod;

    std::mt19937& gen;
    view_t<dtype, 4>& arr;
    const Params& params;
    bool _mapq_mask;

public:

    const seq_t& reference;
    seq_t query;
    TinyHTS::CIGAR cigar;
    TinyHTS::PHRED phred;
    int mapq;
    hts_pos_t pos;

    Record(
        const seq_t& reference,
        const Params& params,
        view_t<dtype, 4>& arr,
        std::mt19937& gen
    ) : gen(gen), arr(arr), params(params), reference(reference) {

        hts_pos_t length = reference.size();

        // Generate a random alignment position
        pos = _random_int(gen, 0, length - 2);

        // Generate random lengths of insertions and deletions
        _ins = _random_int(gen, 0, length);
        _del = _random_int(gen, 0, length - pos - 2);
        _match = length - pos - _del;

        // Start at the 3' end of the sequence
        _rpos = length;
        _qpos = _match + _ins;
        _last_mod = _rpos + params.collapse;

        // Generate a random mapping quality
        mapq = _random_mapq(gen);

        // Pretend the read is low-quality if it is too short or too long
        if (_qpos < params.min_length || _qpos > params.max_length) {
            _mapq_mask = false;
        } else {
            _mapq_mask = (mapq >= params.min_mapq);
        }

        // Generate a random PHRED score per base
        phred = _random_phred(gen, _qpos);

        // Generate the query sequence
        while (_match > 0 || _ins > 0 || _del > 0) { extend(); }

        // Reverse as we generated 3' -> 5'
        std::reverse(query.begin(), query.end());
        std::reverse(cigar.begin(), cigar.end());

    }

    dtype mask(hts_pos_t qpos) const {

        bool _phred_mask = phred.check(qpos, params.min_phred, params.quality_window);
        return static_cast<dtype>(_mapq_mask && _phred_mask);

    }

    void match() {

        _qpos--;
        _rpos--;
        base_t rbase = reference[_rpos];
        query.push_back(rbase);

        arr(_rpos, rbase, rbase) += mask(_qpos);

    }

    void mismatch(bool mod = true) {

        _qpos--;
        _rpos--;
        base_t rbase = reference[_rpos];
        base_t qbase = _random_base_excluding(gen, rbase);
        query.push_back(qbase);

        if (mod && _last_mod - _rpos >= params.collapse && params.mismatches) {
            arr(_rpos, rbase, qbase) += mask(_qpos);
            _last_mod = _rpos;
        } else {
            arr(_rpos, rbase, rbase) += mask(_qpos);
        }

    }

    void ins(bool mod = true) {

        _qpos--;
        base_t qbase = _random_base(gen);
        query.push_back(qbase);

        if (mod && _last_mod - _rpos >= params.collapse && params.insertions) {
            arr(_rpos, qbase, IX_INS) += mask(_qpos);
            _last_mod = _rpos;
        }

    }

    void del(bool mod = true) {

        _rpos--;
        base_t rbase = reference[_rpos];

        if (mod && _last_mod - _rpos >= params.collapse && params.deletions) {
            arr(_rpos, rbase, IX_DEL) += mask(_qpos);
            _last_mod = _rpos;
        } else {
            arr(_rpos, rbase, rbase)  += mask(_qpos);
        }

    }

    void extend() {

        TinyHTS::CIGAR_op op =_random_cigar_op(gen, _match, _ins, _del, cigar.back().type());
        cigar.extend(op);

        switch (op.type()) {

            case TinyHTS::CIGAR_t::MATCH: {

                while (op.advance()) { match(); }

                _match -= op.length();
                break;

            }

            case TinyHTS::CIGAR_t::MISMATCH: {

                mismatch();
                op.advance();
                while (op.advance()) { mismatch(false); }

                _match -= op.length();
                break;

            }

            case TinyHTS::CIGAR_t::DEL: {

                if (op.length() <= params.max_indel_length) {
                    del();
                    op.advance();
                }
                while (op.advance()) { del(false); }

                _del -= op.length();
                break;

            }

            case TinyHTS::CIGAR_t::INS: {

                if (op.length() <= params.max_indel_length && _rpos < reference.size()) {
                    ins();
                    op.advance();
                }
                while (op.advance()) { ins(false); }

                _ins -= op.length();
                break;

            }

            default: {

                throw std::runtime_error("Error generating test case.");

            }

        }

    }

    std::string str(const std::string& ref_name) const {

        int flag = 0;
        std::string rnext = "*";
        int pnext = 0;
        int tlen = 0;
        std::string read_name = "read";

        std::stringstream ss;
        ss << read_name << "\t" << flag << "\t" << ref_name << "\t" << (pos + 1) << "\t" << mapq << "\t" << cigar.str() << "\t" << rnext << "\t" << pnext << "\t" << tlen << "\t" << TinyHTS::_to_str(query) << "\t" << phred.str();

        return ss.str();

    }

};


template <typename dtype>
std::vector<Record<dtype>> _random_records(
    const seq_t& reference,
    const Params& params,
    view_t<dtype, 4>& arr,
    hts_pos_t count,
    std::mt19937& gen
) {

    std::vector<Record<dtype>> records;

    for (hts_pos_t ix = 0; ix < count; ix++) {
        Record<dtype> record(reference, params, arr, gen);
        records.push_back(record);
    }

    return records;

}


static inline std::string __sam_header_sorted() {

    return "@HD\tVN:1.6\tSO:unsorted";

}

static inline std::string __sam_header_ref(size_t n, size_t length) {

    std::stringstream ss;
    ss << "@SQ\tSN:ref" << n << "\tLN:" << length;

    return ss.str();

}

static inline void _write_header(size_t references, size_t length, std::ostream& file) {

    file << __sam_header_sorted() << "\n";
    for (size_t ix = 0; ix < references; ix++) {
        file << __sam_header_ref(ix, length) << "\n";
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
    size_t references,
    size_t queries,
    size_t length,
    Params params,
    std::string out_sam,
    std::string out_fasta,
    std::string out_h5
) {

    std::mt19937 gen = _init_gen();

    std::optional<HDF5::File> _hdf5;
    try {
        _hdf5.emplace(out_h5, HDF5::RWC, mpi, references, 3);
    } catch (const std::exception& e) {
        mpi.err() << "Error: " << e.what() << "\n";
    }
    HDF5::File hdf5 = std::move(_hdf5.value());

    std::ofstream sam(out_sam);
    _write_header(references, length, sam);

    OldSchool::FASTA fasta(out_fasta);

    std::vector<size_t> dims = {references, static_cast<size_t>(length), 4, 6};
    HDF5::Memspace memspace = hdf5.memspace<float, 4>(dims, _path(out_sam));

    std::vector<Record<float>> records;
    for (size_t ix = 0; ix < references; ix++) {

        seq_t reference = _random_sequence(length, gen);
        std::string name = "ref" + std::to_string(ix);
        fasta.write(name, TinyHTS::_to_str(reference));

        view_t<float, 4> arr = memspace.view(ix);
        std::vector<Record<float>> records = _random_records<float>(reference, params, arr, queries, gen);

        for (const auto& record : records) {
            sam << record.str(name) << "\n";
        }

    }

    memspace.safe_write(0);

}

int main(int argc, char** argv) {

    // Initialize MPI processes
    MPI::Manager mpi(argc, argv);
    // Initialize HDF5 manager
    HDF5::Manager _hdf5_manager;
    // Disable native HTS logging as it does not work well in the
    // multi-threaded environment
    TinyHTS::__disable_logging();

    // Parse command line arguments
    testsProgram opt;
    try {
        opt.parse(argc, argv);
    } catch (const std::exception& err) {
        mpi.err() << "Error: " << err.what() << "\n";
        return EXIT_FAILURE;
    }

    Params params = {
        opt.min_mapq,
        opt.min_phred,
        opt.min_length,
        opt.max_length,
        opt.max_indel_length,
        opt.quality_window,
        opt.collapse,
        !opt.no_mismatch,
        !opt.no_insertion,
        !opt.no_deletion
    };

    try {
        __run_tests(mpi, opt.references, opt.queries, opt.length, params, opt.out_sam, opt.out_fasta, opt.out_h5);
    } catch (const std::exception& err) {
        mpi.err() << "Error: " << err.what() << "\n";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;

}
