#ifndef _CMUTS_IO_COMMON_ALIGNMENT_HEADER
#define _CMUTS_IO_COMMON_ALIGNMENT_HEADER

#include <cstdint>
#include <vector>
#include <string>
#include <memory>
#include <iostream>
#include <iomanip>

#include "types.hpp"
#include "stream.hpp"


namespace HTS {


//
// Enums
//


enum class FileType : uint8_t {

    SAM,
    BAM,
    CRAM

};


enum class CIGAR_t : uint8_t {

    // HTS operations

    UNKNOWN,  // Unknown
    MATCH,    // Match
    DEL,      // Deletion
    INS,      // Insertion
    SOFT,     // Soft clipping
    HARD,     // Hard clipping
    SKIP,     // ?
    PAD,      // ?
    BACK,     // ?

    // Additional cmuts operations

    MISMATCH, // Mismatch
    TERM      // Termination

};


//
// Base conversion utilities
//


void _disable_logging();
base_t _from_char(char base);
base_t _from_str(const std::string& base);
std::string str(const seq_t& sequence, int32_t start, int32_t end);
std::string str(const seq_t& sequence);


//
// Header-related utilities
//


struct Header {

    bool sorted        = false;
    int32_t references = 0;

};


FileType _get_filetype(BGZF* _bgzf_file);
Header _read_sam_header(std::unique_ptr<ByteStream>& block, int32_t length, bool read = true);
Header _read_sam_header(std::unique_ptr<ByteStream>& block, bool read = true);


//
// CIGAR operations and MD tags
//


class CIGAR_op {
private:

    CIGAR_t _type   = CIGAR_t::UNKNOWN;
    int32_t _length = 0;
    int32_t _pos    = 0;

public:

    // For susbtitutions and insertions, this will hold the query
    // base at the 3'-most end of the op, which is all that is necessary
    // for cmuts.

    base_t base = IX_UNK;
    std::vector<base_t> sub;

    base_t last() const;
    base_t last(base_t reference) const;

    CIGAR_op() = default;
    explicit CIGAR_op(CIGAR_t type);
    CIGAR_op(CIGAR_t type, int32_t length);

    // Get the CIGAR op
    CIGAR_t type() const;

    // Get the length of the op
    int32_t length() const;

    // Get the length of the query the op consumes
    int32_t qlength() const;

    // Get the length of the reference the op consumes
    int32_t rlength() const;

    // The op as it would appear in a SAM/BAM file
    std::string str() const;

    // Increase the length by val
    void extend(int32_t val);

    // Check if the op has length 0
    bool empty() const;

    // Split the op in two, one with length n and one with length (_length - n)
    std::pair<CIGAR_op, CIGAR_op> split(int32_t n) const;

    // Move forward by 1 in the op, or return false if at the end. Finish -offset early,
    // so e.g. advance(-1) advances to the second-to-last base.
    bool advance(int32_t offset = 0);

    // Return a MATCH of the same length
    CIGAR_op match() const;

};


class CIGAR {
private:

    // The ops in the CIGAR
    std::vector<CIGAR_op> _str;

public:

    // Create an empty CIGAR; possibly with reserved memory
    CIGAR() = default;
    explicit CIGAR(int32_t size);

    // The number of ops

    int32_t size() const;

    // The total number of bases covered by the CIGAR
    // (matches + mismatches + insertions + deletions)

    int32_t bases() const;

    // Get the length of the reference the string consumes

    int32_t rlength() const;

    // Get the length of the query the string consumes

    int32_t qlength() const;

    // Add an op to the end of the CIGAR.
    // Merge it into the final op if it is of the same type.

    void extend(const CIGAR_op& op);

    // Add an op to the end of the CIGAR

    void append(const CIGAR_op& op);

    // Get the final op, or UNKNOWN if empty

    CIGAR_op back() const;

    // The CIGAR as it would appear in a SAM file

    std::string str() const;

    // Is the CIGAR empty?

    bool empty() const;

    // The Hamming distance to the reference sequence

    int32_t hamming() const;

    // Collapse neighbouring modificaitons

    CIGAR revcol(int32_t window) const;

    // Iterator stuff

    auto begin() -> decltype(_str.begin());
    auto begin() const -> decltype(_str.begin());
    auto end() -> decltype(_str.end());
    auto end() const -> decltype(_str.end());

};


//
// PHRED
//


template <typename dtype>
class BaseMask {
public:

    std::vector<dtype> mask;
    int32_t good = 0;

    BaseMask(int32_t length);

};


class PHRED {
private:

    std::vector<qual_t> _qualities;

public:

    PHRED() = default;
    explicit PHRED(const std::vector<qual_t>& qualities);

    qual_t operator[](int32_t ix) const;
    std::string str() const;
    bool check(int32_t ix, qual_t min, int32_t window) const;

    template <typename dtype>
    BaseMask<dtype> mask(qual_t min, int32_t window) const;

};


//
// Alignment
//


class Alignment {
public:

    bool aligned  = false;
    bool primary  = false;
    bool reversed = false;

    qual_t mapq       = MISSING_MAPQ;
    int32_t length    = 0;
    int32_t offset    = 0;
    int32_t reference = -1;

    CIGAR cigar;
    PHRED phred;

    /// Default constructor for empty/unaligned alignments
    Alignment() = default;

    /// Named factory method - makes field order explicit and prevents
    /// aggregate initialization bugs (field order is self-documenting)
    static Alignment create(
        bool aligned,
        bool primary,
        bool reversed,
        qual_t mapq,
        int32_t length,
        int32_t offset,
        int32_t reference,
        CIGAR cigar,
        PHRED phred
    ) {
        Alignment aln;
        aln.aligned   = aligned;
        aln.primary   = primary;
        aln.reversed  = reversed;
        aln.mapq      = mapq;
        aln.length    = length;
        aln.offset    = offset;
        aln.reference = reference;
        aln.cigar     = std::move(cigar);
        aln.phred     = std::move(phred);
        return aln;
    }

};


//
// Printing overloads
//


std::ostream& operator<<(std::ostream& os, FileType fileType);
std::ostream& operator<<(std::ostream& os, CIGAR_t cigar);
std::ostream& operator<<(std::ostream& os, const CIGAR_op& op);
std::ostream& operator<<(std::ostream& os, const CIGAR& cigar);


} // namespace HTS


#endif
