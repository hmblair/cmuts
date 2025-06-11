#ifndef _CMUTS_COMMON_HEADER
#define _CMUTS_COMMON_HEADER

#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <cstring>
#include <sys/types.h>
#include <vector>
#include <fstream>
#include <unordered_set>
#include <span>
#include <filesystem>
#include <algorithm>
extern "C" {
    #include <htslib/sam.h>
    #include <htslib/hts.h>
    #include <htslib/hts_log.h>
    #include <htslib/bgzf.h>
    #include <htslib/kstring.h>
    #include <zlib.h>
    #include <rANS_static.h>
}

// Canonical bases

const int32_t BASES = 4;

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

// BGZF/header constants

const int32_t     BGZF_BUFFER     = 1L << 16;
const std::string BGZF_SORTED     = "coordinate";
const std::string HEADER_SORT_KEY = "SO:";
const std::string HEADER_LEN_KEY  = "LN:";
const int32_t MAGIC_SIZE = 4;
const std::string SAM_MAGIC  = "SAM\1";
const std::string BAM_MAGIC  = "BAM\1";
const std::string CRAM_MAGIC = "CRAM";
const std::string CMUTS_INDEX = ".cmix";

const int32_t ZLIB_WINDOW = 15 + 32;
const int32_t UNALIGNED = -1;

// Quality score special values

const qual_t MAX_MAPQ     = 254;
const qual_t MISSING_MAPQ = 255;
const qual_t MAX_PHRED    = 41;

// Conversion to PHRED ASCII

const int32_t PHRED_OFFSET = 33;

// Reserve this many CIGAR ops when creating a new CIGAR

const int32_t CIGAR_RESERVE = 64;
const int32_t LINE_RESERVE  = 512;

// ITF8 constants

const int32_t ITF8_MAX_SIZE = 5;
const std::array<uint8_t, 5> ITF8_MASK = {0x7F, 0x3F, 0x1F, 0x0F, 0x0F};

// Bit twiddling stuff

constexpr int32_t BYTE = CHAR_BIT;

constexpr std::array<uint8_t, BYTE + 1> _bit_mask() {

    constexpr int32_t _MASK = (1 << BYTE) - 1;
    std::array<uint8_t, BYTE + 1> mask{};

    for (int32_t ix = 0; ix <= BYTE; ix++) {
        mask[ix] = (_MASK << (BYTE - ix)) & _MASK;
    }

    return mask;

}

const std::array<uint8_t, BYTE + 1> BIT_MASK = _bit_mask();




//
// General file utilities
//





bool _exists(const std::string& filename);
void _throw_if_exists(const std::string& filename);
void _throw_if_not_exists(const std::string& filename);
void _safe_move(const std::string& src, const std::string& dst);
std::string _stem(const std::string& name);
std::string _path(const std::string& name);
bool _has_duplicate_paths(const std::vector<std::string>& paths);
void _throw_if_has_duplicate_paths(const std::vector<std::string>& paths);





//
// BGZF file utilities
//





BGZF* _open_bgzf(const std::string& name);
void _close_bgzf(BGZF*& _bgzf_file);
void _read_bgzf(BGZF* _bgzf_file, void* buffer, int32_t size);
std::vector<unsigned char> _read_bgzf(BGZF* _bgzf_file, int32_t size);
std::string _read_bgzf_line(BGZF* _bgzf_file, char* buffer, int32_t& ix, int32_t length);
void _seek_bgzf(BGZF* _bgzf_file, int64_t ptr);
template <typename dtype>
dtype _read_bgzf_single(BGZF* _bgzf_file);





//
// FASTA
//





class FASTA {
private:

    std::string _name;
    std::fstream _file;

public:

    explicit FASTA(const std::string& name);
    std::string next();
    void write(const std::string& name, const std::string& sequence);

};





//
// ByteStream
//





class ByteStream {
protected:

    int32_t _size      = 0;
    int32_t _remaining = 0;

    uint8_t _byte   = 0;
    uint8_t _n_bits = 0;

public:

    explicit ByteStream(int32_t size);
    virtual ~ByteStream() = default;
    ByteStream(const ByteStream&) = delete;
    ByteStream& operator=(const ByteStream&) = delete;

    virtual uint8_t byte()                               = 0;
    virtual std::vector<uint8_t> bytes(int32_t length)   = 0;
    virtual std::vector<uint8_t> line(uint8_t delimiter) = 0;
    virtual std::string str(uint8_t delimiter = '\n')    = 0;

    int32_t size() const;
    bool end() const noexcept;

    uint8_t bits(int32_t length);
    void memcpy(uint8_t* dest, int32_t n);
    int32_t itf8();
    int64_t ltf8();

};





//
// Streams from in-memory data
//





class dataStream : public ByteStream {
private:

    std::span<const uint8_t> _data;
    int32_t _pos = 0;

public:

    explicit dataStream(std::span<const uint8_t> data);

    uint8_t byte() override;
    std::vector<uint8_t> bytes(int32_t length) override;
    std::vector<uint8_t> line(uint8_t delimiter) override;
    std::string str(uint8_t delimiter) override;


};


class ransStream : public dataStream {
private:

    std::vector<uint8_t> _rans_data;
    explicit ransStream(std::vector<uint8_t> data);

public:

    explicit ransStream(std::span<const uint8_t> data, int32_t raw);

};


class zlibStream : public ByteStream {
protected:

    z_stream _stream;
    std::vector<uint8_t> _buffer;

    int32_t _buffer_pos = 0;
    int32_t _buffer_end = 0;
    bool _open          = false;

public:

    zlibStream(std::span<const uint8_t> data, int32_t raw, int32_t buffer);
    ~zlibStream() override;

    virtual void fill();

    uint8_t byte() override;
    std::vector<uint8_t> bytes(int32_t length) override;
    std::vector<uint8_t> line(uint8_t delimiter) override;
    std::string str(uint8_t delimiter) override;

    void update(std::span<const uint8_t> data);

};






//
// Streams from in-storage data
//





class bgzfFileStream : public ByteStream {
private:

    BGZF* _file;

public:

    bgzfFileStream(BGZF* file, int32_t size);

    uint8_t byte() override;
    std::vector<uint8_t> bytes(int32_t length) override;
    std::vector<uint8_t> line(uint8_t delimiter) override;
    std::string str(uint8_t delimiter) override;

};


class zlibFileStream : public zlibStream {
private:

    bgzfFileStream _bgzf;
    std::vector<uint8_t> _bgzf_buffer;
    std::span<const uint8_t> _buffer_span;

    int32_t _bgzf_buffer_pos = 0;
    int32_t _in_remaining    = 0;

public:

    zlibFileStream(BGZF* file, int32_t size, int32_t raw, int32_t buffer);

    void fill() override;

};






namespace HTS {





void _disable_logging();
base_t _from_char(char base);
std::string str(const seq_t& sequence);





//
// Enums
//





enum class FileType : uint8_t {

    SAM,
    BAM,
    CRAM

};


enum class CIGAR_t : uint8_t {

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
// Header-related utilities
//



struct Header {

    bool sorted        = false;
    int32_t references = 0;

};


FileType _get_filetype(BGZF* _bgzf_file);
Header _read_sam_header(std::unique_ptr<ByteStream>& block);






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

    // Get the length of the reference the op consumes

    int32_t rlength() const;

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
    explicit PHRED(const std::vector<qual_t>& qualities);

    qual_t operator[](int32_t ix) const;
    std::string str() const;
    bool check(int32_t ix, qual_t min, int32_t window) const;

    template <typename dtype>
    std::vector<dtype> mask(qual_t min, int32_t window) const;

};





//
// Alignment
//





class Alignment {
public:

    bool aligned   = false;
    qual_t mapq    = MISSING_MAPQ;
    int32_t length = 0;
    int32_t offset = 0;
    int32_t reference = -1;
    CIGAR cigar;
    PHRED phred;

};



class Iterator {
protected:

    int64_t _reads = 0;
    int64_t _curr  = 0;

public:

    explicit Iterator(int64_t reads);
    virtual ~Iterator() = default;
    Iterator(const Iterator&) = delete;
    Iterator& operator=(const Iterator&) = delete;
    Iterator(Iterator&&) = delete;
    Iterator& operator=(Iterator&&) = delete;

    virtual Alignment next() = 0;
    bool end() const noexcept;

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
    int32_t _references;

public:

    Index() = default;
    Index(const std::string& filename, int32_t references);

    IndexBlock read(int32_t ix);
    int64_t aligned();
    int64_t unaligned();

};


Index _sam_index(BGZF* file);
Index _bam_index(BGZF* file);
Index _cram_index(BGZF* file);





//
// File
//





class File {
protected:

    std::string _name;
    FileType _type;

    BGZF* _hts_bgzf     = nullptr;
    bool _sorted        = false;
    int32_t _references = 0;

    Index _index;

public:

    explicit File(const std::string& name, FileType type);
    File(File&& other) noexcept;
    File& operator=(File&& other) = delete;
    File(const File&) = delete;
    File& operator=(const File&) = delete;
    virtual ~File();

    std::string name() const;
    FileType type() const;
    int32_t size() const;
    int64_t aligned();
    int64_t unaligned();

    virtual std::shared_ptr<Iterator> get(int32_t ix, bool seek) = 0;

};


std::unique_ptr<File> _get_sam(const std::string& name);
std::unique_ptr<File> _get_bam(const std::string& name);
std::unique_ptr<File> _get_cram(const std::string& name);





//
// FileGroup
//





class FileGroup {
private:

    std::vector<std::unique_ptr<File>> _group;

public:

    explicit FileGroup(const std::vector<std::string>& filenames);

    int32_t size() const;
    int64_t aligned();
    int64_t unaligned();
    int32_t references() const;

    auto begin() -> decltype(_group.begin());
    auto begin() const -> decltype(_group.begin());
    auto end() -> decltype(_group.end());
    auto end() const -> decltype(_group.end());

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
