#ifndef HTS_HEADER
#define HTS_HEADER

#include <stdexcept>
#include <string>
#include <filesystem>
#include <unordered_set>
#include <span>
#include <fstream>
#include <iostream>
#include <random>
extern "C" {
    #include <htslib/sam.h>
    #include <htslib/hts.h>
    #include <htslib/faidx.h>
    #include <htslib/hts_log.h>
    #include <htslib/kstring.h>
}
#include "mpi.hpp"
#include "hdf5.hpp"

// Canonical bases
const char A = 'A';
const char C = 'C';
const char G = 'G';
const char T = 'T';
const char U = 'U';
const char N = 'N';
// Indices of the bases in HTS files
const uint8_t HTS_A = 1;
const uint8_t HTS_C = 2;
const uint8_t HTS_G = 4;
const uint8_t HTS_T = 8;
const uint8_t HTS_N = 15;
// Internal types used to represent bases and sequences
using base_t = int8_t;
using seq_t  = std::vector<base_t>;
// Indices of the bases in the output arrays
const base_t IX_A   = 0;
const base_t IX_C   = 1;
const base_t IX_G   = 2;
const base_t IX_U   = 3;
const base_t IX_T   = 3;
const base_t IX_DEL = 4;
const base_t IX_INS = 5;
const base_t IX_UNK = -1;
// Index file extensions
const std::string FASTA_INDEX = ".fai";
const std::string BAM_INDEX   = ".bai";
const std::string CRAM_INDEX  = ".crai";
// MD tag constants
const char MD_DEL  = '^';
const char MD_NULL = '0';
// Conversion to PHRED ASCII
const hts_pos_t PHRED_OFFSET = 33;
// Quality score special values
const uint8_t MAX_MAPQ     = 254;
const uint8_t MISSING_MAPQ = 255;
const uint8_t MAX_PHRED    = 41;
// HDF5 dataset to save FASTA sequences in
const std::string SEQUENCE_DS = "sequence";
const size_t DS_DIMS = 2;

namespace HTS {



//
// Enums
//



enum class FileType {
    SAM,
    BAM,
    CRAM
};
std::ostream& operator<<(std::ostream& os, FileType fileType);

enum class CIGAR_t {
    UNKNOWN,  // Unknown
    MATCH,    // Match
    MISMATCH, // Mismatch
    DEL,      // Deletion
    INS,      // Insertion
    SOFT,     // Soft clipping
    HARD,     // Hard clipping
    SKIP,     // ?
    PAD,      // ?
    BACK      // ?
};
std::ostream& operator<<(std::ostream& os, CIGAR_t cigar);



//
// Miscellaneous functions
//



void __disable_logging();

std::string _to_str(const seq_t& sequence);

void _default_write_to_fasta(
    const std::string& filename,
    const std::vector<std::string>& sequences
);



//
// PHRED
//



class PHRED {
private:

    std::vector<uint8_t> _qualities;

public:

    PHRED() = default;
    PHRED(const std::vector<uint8_t>& qualities);
    uint8_t operator[](size_t ix) const;
    std::string str() const;
    bool check(size_t ix, uint8_t min, size_t window) const;

};



//
// FASTA
//



class FASTA {
private:

    std::string _name;
    faidx_t* _hts_fai = nullptr;
    hts_pos_t _size   = 0;

public:

    explicit FASTA(const std::string& filename);
    FASTA(FASTA&& other) noexcept;
    FASTA& operator=(FASTA&& other) noexcept;
    FASTA(const FASTA&);
    FASTA& operator=(const FASTA&) = delete;
    ~FASTA();
    // Get a pointer to the fai object
    faidx_t* ptr() const;
    // Get the filename
    std::string name() const;
    // Get the number of sequences
    size_t size() const;
    // Get the length of the sequence at position ix
    hts_pos_t length(hts_pos_t ix) const;
    // Get the sequence at position ix
    std::string sequence(hts_pos_t ix) const;
    // Get the sequence at position ix as a vector of integers
    seq_t operator[](hts_pos_t ix) const;
    // Get the sequence at position ix as a vector of integers, inserted into the buffer
    void as_int(hts_pos_t ix, view_t<base_t, DS_DIMS> buffer) const;
    // Get the length of the longest sequence in the file
    hts_pos_t max_length() const;
    // Save in the dataset SEQUENCE_DS inside the given HDF5 file
    void to_hdf5(HDF5::File& hdf5, const MPI::Manager& mpi) const;


};

class CIGAR_op {
private:

    CIGAR_t _type     = CIGAR_t::UNKNOWN;
    hts_pos_t _length = 0;
    hts_pos_t _pos    = 0;

public:

    CIGAR_op() = default;
    CIGAR_op(CIGAR_t type);
    CIGAR_op(CIGAR_t type, hts_pos_t length);
    // Get the CIGAR op
    CIGAR_t type() const;
    // Get the length of the op
    hts_pos_t length() const;
    // Get the length of the reference the op consumes
    hts_pos_t rlength() const;
    // Get the length of the query the op consumes
    hts_pos_t qlength() const;
    // The op as it would appear in an HTS file
    std::string str() const;
    // Increase the length by val
    void extend(hts_pos_t val);
    // Get the current position in the op
    hts_pos_t pos() const;
    // Check if the op has length 0
    bool empty() const;
    // Split the op in two, one with length n and one with length (_length - n)
    std::pair<CIGAR_op, CIGAR_op> split(hts_pos_t n) const;
    // Move forward by 1 in the op, or return false if at the end. Finish -offset early,
    // so e.g. advance(-1) advances to the second-to-last base.
    bool advance(ssize_t offset = 0);
    // Return a MATCH of the same length
    CIGAR_op match() const;

};

class MD {
private:

    // The MD tag
    std::string _tag;
    // The current position in the MD tag
    hts_pos_t pos = 0;
    // The current op
    CIGAR_op curr;
    // Update curr to the op of the current position
    void update();

public:

    // Create an MD object from the raw MD tag
    MD(const std::string& tag);
    // Return the next operation in the MD tag, truncated to a length of max. Store the remainder of the truncation in curr.
    CIGAR_op advance(hts_pos_t max = INT_MAX);
    // Get the raw MD tag
    std::string tag() const;

};

class CIGAR {
private:

    // The ops in the CIGAR
    std::vector<CIGAR_op> _str;

public:

    // Create an empty CIGAR; possibly with reserved memory
    CIGAR() = default;
    CIGAR(size_t size);
    // The number of ops
    size_t size() const;
    // The total number of bases covered by the CIGAR
    // (matches + mismatches + insertions + deletions)
    size_t bases() const;
    // Get the length of the reference the op consumes
    hts_pos_t rlength() const;
    // Add an op to the end of the CIGAR.
    // Merge it into the final op if it is of the same type.
    void extend(CIGAR_op op);
    // Add an op to the end of the CIGAR
    void append(CIGAR_op cigar);
    // Get the op at position ix
    CIGAR_op operator[](size_t ix) const;
    // The CIGAR as it would appear in an HTS file
    std::string str() const;
    // Iterator stuff
    auto begin() -> decltype(_str.begin());
    auto begin() const -> decltype(_str.begin());
    auto end() -> decltype(_str.end());
    auto end() const -> decltype(_str.end());

};
std::ostream& operator<<(std::ostream& os, CIGAR cigar);


class Alignment {
private:

    htsFile*   _hts_file = nullptr;
    bam1_t*    _hts_aln  = nullptr;
    hts_itr_t* _hts_iter = nullptr;

    seq_t _reference;

public:

    Alignment(
        htsFile* _hts_file,
        hts_itr_t* _hts_iter,
        const seq_t& _reference
    );
    ~Alignment();
    Alignment(Alignment&& other) noexcept;
    Alignment& operator=(Alignment&& other) noexcept;
    Alignment(const Alignment&) = delete;
    Alignment& operator=(const Alignment&) = delete;

    base_t qbase(hts_pos_t ix) const;
    base_t rbase(hts_pos_t ix) const;
    hts_pos_t qlength() const;
    hts_pos_t rlength() const;
    const seq_t& reference() const;

    bool next() const;
    int64_t reads() const;

    hts_pos_t offset() const;
    uint8_t mapq() const;
    std::span<const uint8_t> phred() const;
    bool empty() const;

    template <typename dtype>
    std::vector<dtype> mask(uint8_t min_quality, hts_pos_t window) const;

    bool aligned() const;
    CIGAR cigar() const;

};


class File {
private:

    htsFile*   _hts_file   = nullptr;
    hts_idx_t* _hts_index  = nullptr;
    bam_hdr_t* _hts_header = nullptr;

    std::string _name;
    FASTA _fasta;
    const MPI::Manager& _mpi;

public:

    File(const std::string& filename, const FASTA& fasta, const MPI::Manager& mpi);
    File(File&& other) noexcept;
    File& operator=(File&& other) = delete;
    File(const File&) = delete;
    File& operator=(const File&) = delete;
    ~File();

    // The MPI object
    const MPI::Manager& mpi() const;
    // Pointers to the various HTS objects
    htsFile* ptr() const;
    hts_idx_t* index() const;
    bam_hdr_t* header() const;

    // Get another file handle to the same file

    htsFile* handle() const;

    const FASTA& fasta() const;


    // Iterator over sequence ix

    hts_itr_t* iter(hts_pos_t ix) const;

    std::string name() const;
    FileType type() const;

    void reset();
    std::string stem() const;
    std::string path() const;

    seq_t reference(hts_pos_t ix) const;
    int64_t size() const;
    int64_t reads(hts_pos_t ix) const;
    int64_t reads() const;
    int64_t unaligned_reads() const;

    hts_pos_t length(int64_t ix) const;
    hts_pos_t max_length() const;

    Alignment alignment(hts_pos_t ix) const;

};


class FileGroup {
private:

    std::vector<File> group;

public:

    FileGroup(
        const std::vector<std::string>& filenames,
        const FASTA& fasta,
        const MPI::Manager& mpi
    );

    int64_t size() const;
    int64_t reads();
    int64_t unaligned_reads() const;
    int64_t references() const;
    int64_t max_length() const;

    auto begin() -> decltype(group.begin());
    auto begin() const -> decltype(group.begin());
    auto end() -> decltype(group.end());
    auto end() const -> decltype(group.end());

};

}

#endif
