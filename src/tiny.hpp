#ifndef TINY_HEADER
#define TINY_HEADER

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <fstream>
#include <unordered_set>
#include <span>
extern "C" {
    #include <htslib/sam.h>
    #include <htslib/hts.h>
    #include <htslib/hts_log.h>
    #include <htslib/bgzf.h>
}
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
// Indices of the bases in binary sequences
const uint8_t BIN_A   = 0x0;
const uint8_t BIN_C   = 0x1;
const uint8_t BIN_G   = 0x2;
const uint8_t BIN_T   = 0x3;
const uint8_t BIN_UNK = 0x4;
// Other binary sequence data constants
const int64_t BASES_PER_BYTE = 4;
// Internal types used to represent bases, sequences, and qualities
using base_t = int8_t;
using seq_t  = std::vector<base_t>;
using qual_t = uint8_t;
// Indices of the bases in the output arrays
const base_t IX_A   = 0;
const base_t IX_C   = 1;
const base_t IX_G   = 2;
const base_t IX_U   = 3;
const base_t IX_T   = 3;
const base_t IX_DEL = 4;
const base_t IX_INS = 5;
const base_t IX_UNK = -1;
// MD tag constants
const char MD_DEL  = '^';
const char MD_NULL = '0';
// BGZF/header constants
const int64_t     BGZF_BUFFER     = 8192;
const std::string BGZF_SORTED     = "coordinate";
const std::string HEADER_SORT_KEY = "SO:";
const std::string HEADER_LEN_KEY  = "LN:";
const int64_t EOH = -1;
const std::string SAM_FT_STR  = "SAM\1";
const std::string BAM_FT_STR  = "BAM\1";
const std::string CRAM_FT_STR = "CRAM";
// Quality score special values
const qual_t MAX_MAPQ     = 254;
const qual_t MISSING_MAPQ = 255;
const qual_t MAX_PHRED    = 41;
// Conversion to PHRED ASCII
const int64_t PHRED_OFFSET = 33;





namespace OldSchool {





class FASTA {
private:

    std::string _name;
    std::fstream _file;

public:

    FASTA(const std::string& name);
    std::string next();
    void write(const std::string& name, const std::string& sequence);

};





} // namespace OldSchool





namespace TinyHTS {





//
// Enums
//





enum class FileType {
    SAM,
    BAM,
    CRAM
};


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





//
// File and miscellaneous utilities
//





void __disable_logging();
std::string _to_str(const seq_t& sequence);





//
// CIGAR operations and MD tags
//





class CIGAR_op {
private:

    CIGAR_t _type   = CIGAR_t::UNKNOWN;
    int64_t _length = 0;
    int64_t _pos    = 0;

public:

    CIGAR_op() = default;
    CIGAR_op(CIGAR_t type);
    CIGAR_op(CIGAR_t type, int64_t length);
    // Get the CIGAR op
    CIGAR_t type() const;
    // Get the length of the op
    int64_t length() const;
    // Get the length of the query the op consumes
    int64_t qlength() const;
    // Get the length of the reference the op consumes
    int64_t rlength() const;
    // The op as it would appear in a SAM/BAM file
    std::string str() const;
    // Increase the length by val
    void extend(int64_t val);
    // Get the current position in the op
    int64_t pos() const;
    // Check if the op has length 0
    bool empty() const;
    // Split the op in two, one with length n and one with length (_length - n)
    std::pair<CIGAR_op, CIGAR_op> split(int64_t n) const;
    // Move forward by 1 in the op, or return false if at the end. Finish -offset early,
    // so e.g. advance(-1) advances to the second-to-last base.
    bool advance(int64_t offset = 0);
    // Return a MATCH of the same length
    CIGAR_op match() const;

};


class MD_tag {
private:

    // The MD tag
    std::string _tag;
    // The current position in the MD tag
    int64_t pos = 0;
    // The current op
    CIGAR_op curr;
    // Update curr to the op of the current position
    void update();

public:

    // Create an MD_tag object from the raw MD tag
    MD_tag(const std::string& tag);
    // Return the next operation in the MD tag, possibly truncated to a length of max. 
    // Store the remainder of the truncation in curr.
    CIGAR_op advance();
    CIGAR_op advance(int64_t max);
    // The MD tag as it would appear in a SAM file
    std::string str() const;

};


class CIGAR {
private:

    // The ops in the CIGAR
    std::vector<CIGAR_op> _str;

public:

    // Create an empty CIGAR; possibly with reserved memory
    CIGAR() = default;
    CIGAR(int64_t size);
    // The number of ops
    int64_t size() const;
    // The total number of bases covered by the CIGAR
    // (matches + mismatches + insertions + deletions)
    int64_t bases() const;
    // Get the length of the reference the op consumes
    int64_t rlength() const;
    // Add an op to the end of the CIGAR.
    // Merge it into the final op if it is of the same type.
    void extend(CIGAR_op op);
    // Add an op to the end of the CIGAR
    void append(CIGAR_op cigar);
    // Get the op at position ix
    CIGAR_op operator[](int64_t ix) const;
    // Get the final op, or UNKNOWN if empty
    CIGAR_op back() const;
    // The CIGAR as it would appear in a SAM file
    std::string str() const;
    // Iterator stuff
    auto begin() -> decltype(_str.begin());
    auto begin() const -> decltype(_str.begin());
    auto end() -> decltype(_str.end());
    auto end() const -> decltype(_str.end());

};






//
// PHRED
//





class PHRED {
private:

    std::vector<qual_t> _qualities;

public:

    PHRED() = default;
    PHRED(const std::vector<qual_t>& qualities);
    qual_t operator[](int64_t ix) const;
    std::string str() const;
    bool check(int64_t ix, qual_t min, int64_t window) const;

};






//
// Printing overloads
//





std::ostream& operator<<(std::ostream& os, FileType fileType);
std::ostream& operator<<(std::ostream& os, CIGAR_t cigar);
std::ostream& operator<<(std::ostream& os, CIGAR cigar);





//
// FASTA
//





class HeaderBlock{
public:

    int64_t length    = 0;
    int64_t sequences = 0;
    bool empty() const;

};

class FASTA {
private:

    std::string _fasta_name;
    std::string _name;
    std::ifstream _file;
    std::vector<HeaderBlock> _blocks;
    int64_t _offset = 0;

public:

    FASTA(const std::string& name);
    FASTA(FASTA&& other) noexcept;
    FASTA& operator=(FASTA&& other) noexcept;
    FASTA(const FASTA&);
    FASTA& operator=(const FASTA&) = delete;

    std::string name() const;
    int64_t size() const;
    int64_t length(int64_t ix) const;
    seq_t sequence(int64_t ix);
    int64_t longest() const;

};





//
// Index
//





class IndexBlock {
public:

    int64_t ptr   = -1;
    int64_t reads = 0;

    bool empty() const;
    void tell(BGZF* _hts_bgzf);
    void write_ptr(std::ofstream& file);
    void write_reads(std::ofstream& file);

};

class Index {
private:

    std::string _name;
    std::ifstream _file;
    int64_t _references;

public:

    Index() = default;
    Index(const std::string& filename, int64_t references);
    IndexBlock read(int64_t ix);
    int64_t aligned();
    int64_t unaligned();

};





//
// Alignment
//





class Alignment {
private:

    BGZF* _hts_bgzf  = nullptr;
    bam1_t* _hts_aln = nullptr;
    int64_t _reads   = 0;
    int64_t _curr    = 0;

public:

    Alignment(BGZF* _hts_bgzf, bam1_t* _hts_aln, int64_t reads);
    Alignment(Alignment&& other) noexcept;
    Alignment& operator=(Alignment&& other) noexcept;
    Alignment(const Alignment&) = delete;
    Alignment& operator=(const Alignment&) = delete;

    bool next();
    bool aligned() const;
    base_t base(int64_t ix) const;
    int64_t length() const;
    int64_t offset() const;
    qual_t mapq() const;
    std::span<const qual_t> phred() const;
    CIGAR cigar() const;

    template<typename dtype>
    std::vector<dtype> mask(qual_t min, int64_t window) const;

};





//
// File
//





class File {
private:

    std::string _name;
    BGZF* _hts_bgzf;
    bam1_t* _hts_aln;

    FileType _type;
    bool _sorted;
    int64_t _references;
    int64_t _body_ptr;

    Index _index;

public:

    File(const std::string& name);
    File(File&& other) noexcept;
    File& operator=(File&& other) = delete;
    File(const File&) = delete;
    File& operator=(const File&) = delete;
    ~File();

    std::string name() const;
    BGZF* ptr() const;
    FileType type() const;
    int64_t size() const;
    int64_t aligned();
    int64_t unaligned();
    Alignment alignment(int64_t ix, bool seek);

};





//
// FileGroup
//





class FileGroup {
private:

    std::vector<File> _group;

public:

    FileGroup(const std::vector<std::string>& filenames);

    int64_t size() const;
    int64_t aligned();
    int64_t unaligned();
    int64_t references() const;

    auto begin() -> decltype(_group.begin());
    auto begin() const -> decltype(_group.begin());
    auto end() -> decltype(_group.end());
    auto end() const -> decltype(_group.end());

};





} // namespace TinyHTS





#endif
