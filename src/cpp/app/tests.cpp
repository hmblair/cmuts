#include "app/tests.hpp"

namespace TestGen {

//
// Random class implementation
//

Random::Random() {
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    _gen = std::mt19937(seed);
}

Random::Random(unsigned seed) : _gen(seed) {}

int Random::integer(int low, int high) {
    std::uniform_int_distribution<> dist(low, high);
    return dist(_gen);
}

int Random::integer_excluding(int low, int high, int exclude) {
    // Generate in range [low, high-1], then shift if >= exclude
    int value = integer(low, high - 1);
    if (value >= exclude) {
        ++value;
    }
    return value;
}

float Random::uniform() {
    std::uniform_real_distribution<float> dist(0.0f, 1.0f);
    return dist(_gen);
}

bool Random::bernoulli(float probability) {
    return uniform() < probability;
}

base_t Random::base() {
    return static_cast<base_t>(integer(0, BASES - 1));
}

base_t Random::base_excluding(base_t exclude) {
    return static_cast<base_t>(integer_excluding(0, BASES - 1, exclude));
}

qual_t Random::mapq() {
    return static_cast<qual_t>(integer(0, MAX_MAPQ));
}

qual_t Random::phred() {
    return static_cast<qual_t>(integer(0, MAX_PHRED));
}

template <typename T>
T Random::sample(const std::vector<T>& values) {
    int ix = integer(0, static_cast<int>(values.size()) - 1);
    return values[ix];
}

// Explicit instantiation for CIGAR_t
template HTS::CIGAR_t Random::sample(const std::vector<HTS::CIGAR_t>& values);

std::mt19937& Random::generator() {
    return _gen;
}

//
// Sequence generation
//

seq_t random_sequence(int32_t length, Random& rng) {
    seq_t sequence(length);
    for (int32_t ix = 0; ix < length; ix++) {
        sequence[ix] = rng.base();
    }
    return sequence;
}

HTS::PHRED random_qualities(int32_t length, Random& rng) {
    std::vector<qual_t> qualities(length);
    for (auto& q : qualities) {
        q = rng.phred();
    }
    return HTS::PHRED(qualities);
}

//
// AlignmentGenerator implementation
//

AlignmentGenerator::AlignmentGenerator(
    const seq_t& reference,
    const Params& params,
    view_t<float, 4>& expected_output,
    Random& rng
) : _rng(rng),
    _params(params),
    _reference(reference),
    _expected(expected_output)
{
    int32_t ref_length = static_cast<int32_t>(reference.size());

    // Generate random alignment start position
    // Leave room for at least 2 bases of alignment
    _start_pos = rng.integer(0, ref_length - 2);

    // Generate random indel counts
    // Deletions are bounded by remaining reference length
    _ins_remaining = rng.integer(0, ref_length);
    _del_remaining = rng.integer(0, ref_length - _start_pos - 2);
    _match_remaining = ref_length - _start_pos - _del_remaining;

    // Initialize position tracking (building 3'→5')
    // We start at the 3' end and work backwards
    _rpos = ref_length;
    _qpos = _match_remaining + _ins_remaining;
    _last_mod_pos = _rpos + params.collapse;  // Allow first modification

    // Generate mapping quality
    _mapq = rng.mapq();

    // Determine if read passes length filter
    bool passes_length = (_qpos >= params.min_length && _qpos <= params.max_length);

    // Determine if read passes MAPQ filter
    _passes_mapq = passes_length && (_mapq >= params.min_mapq);

    // Randomly make some reads unaligned (for testing unaligned handling)
    _is_aligned = !rng.bernoulli(UNALIGNED_PROBABILITY);
    if (!_is_aligned) {
        _passes_mapq = false;
    }

    // Generate quality scores for the query
    _phred = random_qualities(_qpos, rng);

    // Build the alignment backwards (3'→5')
    // This matches how cmuts processes alignments internally
    while (_match_remaining > 0 || _ins_remaining > 0 || _del_remaining > 0) {
        extend_alignment();
    }

    // Reverse to get correct 5'→3' orientation for output
    std::reverse(_query.begin(), _query.end());
    std::reverse(_cigar.begin(), _cigar.end());
}

float AlignmentGenerator::mask_value(int32_t qpos) const {
    bool passes_phred = _phred.check(qpos, _params.min_phred, _params.quality_window);
    return static_cast<float>(_passes_mapq && passes_phred);
}

HTS::CIGAR_op AlignmentGenerator::select_next_operation() {
    // Handle edge cases where only one operation type remains
    if (_match_remaining == 0) {
        if (_ins_remaining == 0) {
            return HTS::CIGAR_op(HTS::CIGAR_t::DEL, _del_remaining);
        }
        if (_del_remaining == 0) {
            return HTS::CIGAR_op(HTS::CIGAR_t::INS, _ins_remaining);
        }
    }

    // Build list of valid next operations
    // Avoid consecutive operations of the same type
    std::vector<HTS::CIGAR_t> candidates;
    HTS::CIGAR_t prev_type = _cigar.back().type();

    if (_match_remaining > 0) {
        if (prev_type != HTS::CIGAR_t::MATCH) {
            candidates.push_back(HTS::CIGAR_t::MATCH);
        }
        if (prev_type != HTS::CIGAR_t::MISMATCH) {
            candidates.push_back(HTS::CIGAR_t::MISMATCH);
        }
    }
    if (_ins_remaining > 0 && prev_type != HTS::CIGAR_t::INS) {
        candidates.push_back(HTS::CIGAR_t::INS);
    }
    if (_del_remaining > 0 && prev_type != HTS::CIGAR_t::DEL) {
        candidates.push_back(HTS::CIGAR_t::DEL);
    }

    // Select random operation type
    HTS::CIGAR_t type = _rng.sample(candidates);

    // Select random length based on operation type
    int32_t max_length;
    switch (type) {
        case HTS::CIGAR_t::MATCH:
        case HTS::CIGAR_t::MISMATCH:
            max_length = _match_remaining;
            break;
        case HTS::CIGAR_t::INS:
            max_length = _ins_remaining;
            break;
        case HTS::CIGAR_t::DEL:
            max_length = _del_remaining;
            break;
        default:
            CMUTS_THROW("Invalid CIGAR operation type in test generation");
    }

    int32_t length = _rng.integer(1, max_length);
    return HTS::CIGAR_op(type, length);
}

void AlignmentGenerator::add_match() {
    _qpos--;
    _rpos--;
    base_t ref_base = _reference[_rpos];
    _query.push_back(ref_base);

    // Record match in expected output
    _expected(_rpos, ref_base, ref_base) += mask_value(_qpos);
}

void AlignmentGenerator::add_mismatch(bool count_modification) {
    _qpos--;
    _rpos--;
    base_t ref_base = _reference[_rpos];
    base_t query_base = _rng.base_excluding(ref_base);
    _query.push_back(query_base);

    // Check collapse distance and whether mismatches are enabled
    bool should_count = count_modification &&
                        (_last_mod_pos - _rpos >= _params.collapse) &&
                        _params.mismatches;

    if (should_count) {
        _expected(_rpos, ref_base, query_base) += mask_value(_qpos);
        _last_mod_pos = _rpos;
    } else {
        // Count as match if not counting the mismatch
        _expected(_rpos, ref_base, ref_base) += mask_value(_qpos);
    }
}

void AlignmentGenerator::add_insertion(bool count_modification) {
    _qpos--;
    base_t query_base = _rng.base();
    _query.push_back(query_base);

    // Check collapse distance, position validity, and whether insertions are enabled
    bool valid_position = (_rpos <= static_cast<int32_t>(_reference.size()) && _rpos > 0);
    bool should_count = count_modification &&
                        valid_position &&
                        (_last_mod_pos - _rpos >= _params.collapse) &&
                        _params.insertions;

    if (should_count) {
        // Insertions are recorded at the position before the insertion
        _expected(_rpos - 1, query_base, IX_INS) += mask_value(_qpos);
        _last_mod_pos = _rpos;
    }
}

void AlignmentGenerator::add_deletion(bool count_modification) {
    _rpos--;
    base_t ref_base = _reference[_rpos];

    // Check collapse distance and whether deletions are enabled
    bool should_count = count_modification &&
                        (_last_mod_pos - _rpos >= _params.collapse) &&
                        _params.deletions;

    if (should_count) {
        _expected(_rpos, ref_base, IX_DEL) += mask_value(_qpos);
        _last_mod_pos = _rpos;
    } else {
        // Count as match if not counting the deletion
        _expected(_rpos, ref_base, ref_base) += mask_value(_qpos);
    }
}

void AlignmentGenerator::extend_alignment() {
    HTS::CIGAR_op op = select_next_operation();
    _cigar.extend(op);

    switch (op.type()) {
        case HTS::CIGAR_t::MATCH: {
            while (op.advance()) {
                add_match();
            }
            _match_remaining -= op.length();
            break;
        }

        case HTS::CIGAR_t::MISMATCH: {
            // First base of mismatch run counts as modification
            add_mismatch(true);
            op.advance();
            // Subsequent bases don't count as separate modifications
            while (op.advance()) {
                add_mismatch(false);
            }
            _match_remaining -= op.length();
            break;
        }

        case HTS::CIGAR_t::DEL: {
            // Only count if within max indel length
            bool count_first = (op.length() <= _params.max_indel_length);
            if (count_first) {
                add_deletion(true);
                op.advance();
            }
            while (op.advance()) {
                add_deletion(false);
            }
            _del_remaining -= op.length();
            break;
        }

        case HTS::CIGAR_t::INS: {
            // Only count if within max indel length
            bool count_first = (op.length() <= _params.max_indel_length);
            if (count_first) {
                add_insertion(true);
                op.advance();
            }
            while (op.advance()) {
                add_insertion(false);
            }
            _ins_remaining -= op.length();
            break;
        }

        default: {
            CMUTS_THROW("Unexpected CIGAR operation in test generation");
        }
    }
}

const seq_t& AlignmentGenerator::query() const { return _query; }
const HTS::CIGAR& AlignmentGenerator::cigar() const { return _cigar; }
const HTS::PHRED& AlignmentGenerator::phred() const { return _phred; }
int32_t AlignmentGenerator::position() const { return _start_pos; }
qual_t AlignmentGenerator::mapq() const { return _mapq; }
bool AlignmentGenerator::is_aligned() const { return _is_aligned; }

//
// SamRecord implementation
//

SamRecord::SamRecord(
    const std::string& ref_name,
    const AlignmentGenerator& alignment
) : _read_name("read"),
    _flag(0),
    _ref_name(alignment.is_aligned() ? ref_name : "*"),
    _pos(alignment.position() + 1),  // SAM uses 1-based positions
    _mapq(alignment.mapq()),
    _cigar(alignment.cigar().str()),
    _rnext("*"),
    _pnext(0),
    _tlen(0),
    _seq(HTS::str(alignment.query())),
    _qual(alignment.phred().str())
{}

std::string SamRecord::str() const {
    std::stringstream ss;
    ss << _read_name << "\t"
       << _flag << "\t"
       << _ref_name << "\t"
       << _pos << "\t"
       << static_cast<int>(_mapq) << "\t"
       << _cigar << "\t"
       << _rnext << "\t"
       << _pnext << "\t"
       << _tlen << "\t"
       << _seq << "\t"
       << _qual;
    return ss.str();
}

//
// SAM header generation
//

std::string sam_header_line() {
    return "@HD\tVN:1.6\tSO:unsorted";
}

std::string sam_reference_line(size_t index, size_t length) {
    std::stringstream ss;
    ss << "@SQ\tSN:ref" << index << "\tLN:" << length;
    return ss.str();
}

void write_sam_header(size_t num_references, size_t ref_length, std::ostream& out) {
    out << sam_header_line() << "\n";
    for (size_t ix = 0; ix < num_references; ix++) {
        out << sam_reference_line(ix, ref_length) << "\n";
    }
}

//
// Main test generation
//

void generate_test_data(
    const MPI::Manager& mpi,
    size_t num_references,
    size_t reads_per_reference,
    size_t reference_length,
    const Params& params,
    const std::string& sam_path,
    const std::string& fasta_path,
    const std::string& hdf5_path,
    int seed
) {
    // Use provided seed if non-negative, otherwise use time-based seed
    Random rng = (seed >= 0) ? Random(static_cast<unsigned>(seed)) : Random();

    // Determine chunking for HDF5 output
    size_t chunk_size = std::min(DEFAULT_CHUNK_SIZE, num_references);
    size_t num_chunks = (num_references + chunk_size - 1) / chunk_size;

    // Initialize HDF5 output file
    std::optional<HDF5::File> hdf5_opt;
    try {
        hdf5_opt.emplace(hdf5_path, HDF5::RWC, mpi, chunk_size, 3);
    } catch (const std::exception& e) {
        mpi.err() << "Error creating HDF5 file: " << e.what() << "\n";
        return;
    }
    HDF5::File hdf5 = std::move(hdf5_opt.value());

    // Initialize SAM output file
    std::ofstream sam(sam_path);
    write_sam_header(num_references, reference_length, sam);

    // Initialize FASTA output file
    FASTA fasta(fasta_path);

    // Create HDF5 memspace for expected values
    std::vector<size_t> dims = {
        num_references,
        reference_length,
        BASES,
        6  // A, C, G, U, DEL, INS
    };
    HDF5::Memspace memspace = hdf5.memspace<float, 4>(dims, _path(sam_path) + "/counts-1d");

    // Generate test data for each reference
    for (size_t chunk_ix = 0; chunk_ix < num_chunks; chunk_ix++) {
        for (size_t ref_in_chunk = 0; ref_in_chunk < chunk_size; ref_in_chunk++) {
            size_t ref_ix = chunk_ix * chunk_size + ref_in_chunk;
            if (ref_ix >= num_references) break;

            // Generate random reference sequence
            seq_t reference = random_sequence(static_cast<int32_t>(reference_length), rng);
            std::string ref_name = "ref" + std::to_string(ref_ix);

            // Write reference to FASTA
            fasta.write(ref_name, HTS::str(reference));

            // Get view into expected output array for this reference
            view_t<float, 4> expected = memspace.view(ref_in_chunk);

            // Generate alignments for this reference
            for (size_t read_ix = 0; read_ix < reads_per_reference; read_ix++) {
                AlignmentGenerator alignment(reference, params, expected, rng);
                SamRecord record(ref_name, alignment);
                sam << record.str() << "\n";
            }
        }

        // Write chunk to HDF5
        memspace.safe_write(chunk_ix * chunk_size);
        memspace.clear();
    }
}

} // namespace TestGen

//
// Program entry point
//

int main(int argc, char** argv) {
    // Initialize MPI
    MPI::Manager mpi(argc, argv);

    // Initialize HDF5 manager
    HDF5::Manager hdf5_manager;

    // Disable native HTS logging (doesn't work well with multithreading)
    HTS::_disable_logging();

    // Parse command line arguments
    testsProgram opt;
    try {
        opt.parse(argc, argv);
    } catch (const std::exception& err) {
        mpi.err() << "Error: " << err.what() << "\n";
        return EXIT_FAILURE;
    }

    // Build parameters struct from command line options
    TestGen::Params params = {
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

    // Generate test data
    try {
        TestGen::generate_test_data(
            mpi,
            opt.references,
            opt.queries,
            opt.length,
            params,
            opt.out_sam,
            opt.out_fasta,
            opt.out_h5,
            opt.seed
        );
    } catch (const std::exception& err) {
        mpi.err() << "Error: " << err.what() << "\n";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
