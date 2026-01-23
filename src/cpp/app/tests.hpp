#ifndef TESTS_HEADER
#define TESTS_HEADER

#include "core/cmuts.hpp"
#include "io/fasta.hpp"
#include "common.hpp"
#include "infra/utils.hpp"
#include "io/hdf5.hpp"
#include "infra/mpi.hpp"
#include "generated/tests_args.hpp"

#include <random>
#include <chrono>

namespace TestGen {

//
// Constants
//

/// Probability that a read is unaligned (for testing unaligned read handling)
constexpr float UNALIGNED_PROBABILITY = 0.1f;

/// Default HDF5 chunk size for writing test data
constexpr size_t DEFAULT_CHUNK_SIZE = 128;

//
// Parameters for test generation
//

/// Parameters controlling test data generation.
/// These mirror the filtering parameters used by cmuts-core.
struct Params {
    int64_t min_mapq;           ///< Minimum mapping quality threshold
    int64_t min_phred;          ///< Minimum PHRED quality score
    int64_t min_length;         ///< Minimum read length
    int64_t max_length;         ///< Maximum read length
    int64_t max_indel_length;   ///< Maximum indel length to count
    int64_t quality_window;     ///< Window size for quality averaging
    int64_t collapse;           ///< Minimum distance between modifications
    bool mismatches;            ///< Whether to count mismatches
    bool insertions;            ///< Whether to count insertions
    bool deletions;             ///< Whether to count deletions
};

//
// Random number generation utilities
//

/// Encapsulates random number generation for test data.
/// Provides typed random value generation with clear semantics.
class Random {
private:
    std::mt19937 _gen;

public:
    /// Initialize with time-based seed
    Random();

    /// Initialize with specific seed for reproducibility
    explicit Random(unsigned seed);

    /// Generate random integer in [low, high] inclusive
    int integer(int low, int high);

    /// Generate random integer in [low, high] excluding a specific value
    int integer_excluding(int low, int high, int exclude);

    /// Generate random float in [0, 1)
    float uniform();

    /// Generate random boolean with given probability of true
    bool bernoulli(float probability);

    /// Generate random base index (0-3 for A, C, G, T/U)
    base_t base();

    /// Generate random base excluding a specific base
    base_t base_excluding(base_t exclude);

    /// Generate random MAPQ value in valid range
    qual_t mapq();

    /// Generate random PHRED quality score
    qual_t phred();

    /// Sample random element from vector
    template <typename T>
    T sample(const std::vector<T>& values);

    /// Get reference to underlying generator (for compatibility)
    std::mt19937& generator();
};

//
// Sequence generation
//

/// Generate a random nucleotide sequence of given length
seq_t random_sequence(int32_t length, Random& rng);

/// Generate random PHRED quality scores for a read
HTS::PHRED random_qualities(int32_t length, Random& rng);

//
// Alignment generation
//

/// Generates a random alignment and tracks expected output values.
///
/// Alignments are built from 3' to 5' (backwards) to match how cmuts
/// processes alignments internally. The CIGAR string and query sequence
/// are reversed before output.
class AlignmentGenerator {
private:
    // Generation state
    int32_t _match_remaining;    ///< Remaining match/mismatch bases to generate
    int32_t _ins_remaining;      ///< Remaining insertion bases to generate
    int32_t _del_remaining;      ///< Remaining deletion bases to generate
    int32_t _rpos;               ///< Current reference position (decreasing)
    int32_t _qpos;               ///< Current query position (decreasing)
    int32_t _last_mod_pos;       ///< Position of last modification (for collapse)

    // Configuration
    Random& _rng;
    const Params& _params;
    const seq_t& _reference;

    // Alignment properties
    bool _is_aligned;
    bool _passes_mapq;
    int32_t _start_pos;
    qual_t _mapq;

    // Output tracking
    view_t<float, 4>& _expected;
    seq_t _query;
    HTS::CIGAR _cigar;
    HTS::PHRED _phred;

    /// Calculate mask value for current position based on quality filters
    float mask_value(int32_t qpos) const;

    /// Extend alignment with a random CIGAR operation
    void extend_alignment();

    /// Add a match operation
    void add_match();

    /// Add a mismatch operation
    void add_mismatch(bool count_modification);

    /// Add an insertion operation
    void add_insertion(bool count_modification);

    /// Add a deletion operation
    void add_deletion(bool count_modification);

    /// Select next CIGAR operation type based on remaining counts
    HTS::CIGAR_op select_next_operation();

public:
    /// Generate a random alignment against the reference sequence.
    /// Updates the expected output array with counted modifications.
    AlignmentGenerator(
        const seq_t& reference,
        const Params& params,
        view_t<float, 4>& expected_output,
        Random& rng
    );

    /// Get the generated query sequence
    const seq_t& query() const;

    /// Get the generated CIGAR string
    const HTS::CIGAR& cigar() const;

    /// Get the generated PHRED scores
    const HTS::PHRED& phred() const;

    /// Get the alignment start position (0-based)
    int32_t position() const;

    /// Get the mapping quality
    qual_t mapq() const;

    /// Check if the read is aligned
    bool is_aligned() const;
};

//
// SAM output formatting
//

/// Formats alignment data as a SAM record string.
class SamRecord {
private:
    std::string _read_name;
    int _flag;
    std::string _ref_name;
    int32_t _pos;           ///< 1-based position for SAM
    qual_t _mapq;
    std::string _cigar;
    std::string _rnext;
    int _pnext;
    int _tlen;
    std::string _seq;
    std::string _qual;

public:
    /// Construct SAM record from alignment data
    SamRecord(
        const std::string& ref_name,
        const AlignmentGenerator& alignment
    );

    /// Format as SAM record string
    std::string str() const;
};

//
// SAM file header generation
//

/// Generate SAM header line
std::string sam_header_line();

/// Generate SAM reference sequence line
std::string sam_reference_line(size_t index, size_t length);

/// Write complete SAM header to stream
void write_sam_header(size_t num_references, size_t ref_length, std::ostream& out);

//
// Main test generation
//

/// Generate test data files (SAM, FASTA, HDF5 with expected values)
/// @param seed Random seed for reproducibility. Use -1 for time-based seed.
void generate_test_data(
    const MPI::Manager& mpi,
    size_t num_references,
    size_t reads_per_reference,
    size_t reference_length,
    const Params& params,
    const std::string& sam_path,
    const std::string& fasta_path,
    const std::string& hdf5_path,
    int seed = -1
);

} // namespace TestGen

//
// Program entry point
//

class testsProgram : public Program {
public:
    TESTSPROGRAM_ARG_MEMBERS

    testsProgram()
        : Program(TESTSPROGRAM_PROGRAM_NAME, TESTSPROGRAM_PROGRAM_VERSION),
          TESTSPROGRAM_ARG_INIT
    {}
};

#endif
