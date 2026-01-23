#ifndef _CMUTS_IO_COMMON_TYPES_HEADER
#define _CMUTS_IO_COMMON_TYPES_HEADER

#include <cstdint>
#include <string>
#include <vector>
#include <array>
#include <climits>
#include <fstream>

extern "C" {
    #include <htslib/bgzf.h>
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

const base_t IX_A    = 0;
const base_t IX_C    = 1;
const base_t IX_G    = 2;
const base_t IX_U    = 3;
const base_t IX_T    = 3;
const base_t IX_DEL  = 4;
const base_t IX_INS  = 5;
const base_t IX_TERM = 6;
const base_t IX_UNK  = -1;

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

    constexpr int32_t BASE = (1 << BYTE) - 1;
    std::array<uint8_t, BYTE + 1> mask{};

    for (int32_t ix = 0; ix <= BYTE; ix++) {
        mask[ix] = (BASE << (BYTE - ix));
    }

    return mask;

}

const std::array<uint8_t, BYTE + 1> BIT_MASK = _bit_mask();


//
// General file utilities
//


bool _exists(const std::string& filename);
void _delete(const std::string& filename);
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


#endif
