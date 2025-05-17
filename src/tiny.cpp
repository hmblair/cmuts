#include "tiny.hpp"
#include <fstream>





static inline bool _exists(const std::string& filename) {

    return std::filesystem::exists(filename);

}


static inline void _throw_if_exists(const std::string& filename) {

    if (_exists(filename)) {
        throw std::runtime_error("The file \"" + filename + "\" already exists.");
    }

}


static inline void _throw_if_not_exists(const std::string& filename) {

    if (!_exists(filename)) {
        throw std::runtime_error("The file \"" + filename + "\" does not exist.");
    }

}


static inline void _safe_move(const std::string& src, const std::string& dst) {

    _throw_if_not_exists(src);
    _throw_if_exists(dst);
    std::filesystem::rename(src, dst);

}


static inline std::string _stem(const std::string& name) {

    std::filesystem::path path(name);
    return path.stem().string();

}


static inline std::string _path(const std::string& name) {

    std::filesystem::path path(name);
    std::string stem = path.stem().string();
    std::string parent = path.parent_path().string();
    if (parent.empty()) { return stem; }
    return parent + "/" + stem;

}


static inline bool _has_duplicate_paths(const std::vector<std::string>& paths) {

    std::unordered_set<std::string> seen;
    for (const auto& path : paths) {
        std::string transformedPath = _path(path);
        if (seen.find(transformedPath) != seen.end()) {
            return true;
        }
        seen.insert(transformedPath);
    }
    return false;

}


static inline void _throw_if_has_duplicate_paths(const std::vector<std::string>& paths) {

    if (_has_duplicate_paths(paths)) {
        throw std::runtime_error("Duplicate input paths (modulo extensions) are not allowed.");
    }

}





namespace OldSchool {





//
// Reading and writing regular FASTA files
//





static inline bool _is_fasta_header(const std::string& line) {

    return !line.empty() && line[0] == '>';

}


static inline std::string _get_fasta_sequence(std::fstream& file) {

    std::string line, sequence;
    while (std::getline(file, line) && !_is_fasta_header(line)) {
        sequence += line;
    }

    return sequence;

}


static inline void _write_fasta_record(const std::string& name, const std::string& sequence, std::fstream& file) {

    file << ">" << name << "\n";
    file << sequence << "\n";

}


FASTA::FASTA(const std::string& name) : _name(name) {

    if (!_exists(name)) {
        _file.open(name, std::ios::out);
    } else {
        _file.open(name, std::ios::in);
        std::string line;
        std::getline(_file, line);
    }

}


std::string FASTA::next() {

    return _get_fasta_sequence(_file);

}


void FASTA::write(const std::string& name, const std::string& sequence) {

    _write_fasta_record(name, sequence, _file);

}





} // namespace OldSchool





namespace TinyHTS {





//
// File and miscellaneous utilities
//





void __disable_logging() { hts_set_log_level(htsLogLevel::HTS_LOG_OFF); }





//
// BGZF tools
//





static inline BGZF* _open_bgzf(const std::string& name) {

    BGZF* _bgzf_file = bgzf_open(name.c_str(), "r");
    if (_bgzf_file == nullptr) {
        throw std::runtime_error("Failed to open the file \"" + name + "\" using BGZF.");
    }

    return _bgzf_file;

}


static inline void _close_bgzf(BGZF*& _bgzf_file) {

    if (_bgzf_file != nullptr) {
        bgzf_close(_bgzf_file);
        _bgzf_file = nullptr;
    }

}


static inline void _read_bgzf(BGZF* _bgzf_file, void* buffer, int64_t size) {

    int64_t bytes = bgzf_read(_bgzf_file, buffer, size);
    if (bytes != size) {
        _close_bgzf(_bgzf_file);
        throw std::runtime_error("BGZF failed to read " + std::to_string(size - bytes) + " of " + std::to_string(size) + " bytes.");
    }

}


template <typename dtype>
static inline dtype _read_bgzf_single(BGZF* _bgzf_file) {

    dtype value;
    _read_bgzf(_bgzf_file, &value, sizeof(dtype));
    return value;

}


static inline FileType _get_filetype(BGZF* _bgzf_file) {

    char magic[4];
    _read_bgzf(_bgzf_file, magic, 4);

    if (strncmp(magic, SAM_FT_STR.c_str(), 4) == 0) {

        return FileType::SAM;

    } else if (strncmp(magic, BAM_FT_STR.c_str(), 4) == 0) {

        return FileType::BAM;

    } else if (strncmp(magic, CRAM_FT_STR.c_str(), 4) == 0) {

        return FileType::CRAM;

    } else {

        _close_bgzf(_bgzf_file);
        throw std::runtime_error("Magic check failed -- invalid BAM file.");

    }

}


static inline void _seek_bgzf(BGZF* _bgzf_file, int64_t ptr) {

    if (bgzf_seek(_bgzf_file, ptr, SEEK_SET) < 0) {

        _close_bgzf(_bgzf_file);
        throw std::runtime_error("Failed to seek to position " + std::to_string(ptr) + ".");

    }

}


static inline std::string _read_bgzf_line(
    BGZF* _bgzf_file,
    char* buffer,
    int64_t& ix,
    int64_t length
) {

    int len = 0;
    while (len < BGZF_BUFFER - 1 && ix < length) {

        if (bgzf_read(_bgzf_file, &buffer[len], 1) != 1) { break; }
        ix++;

        if (buffer[len] == '\n') { break; }
        len++;

    }

    return std::string(buffer, len);

}






//
// Converting to and from internal HTSlib representaitons
//







static inline base_t _base_from_hts(uint8_t base) {

    switch (base) {
        case HTS_A: { return IX_A;   }
        case HTS_C: { return IX_C;   }
        case HTS_G: { return IX_G;   }
        case HTS_T: { return IX_T;   }
        default:    { return IX_UNK; }
    }

}


static inline CIGAR_op _cigar_from_hts(uint32_t op) {

    int64_t length = bam_cigar_oplen(op);
    CIGAR_t type;

    switch (bam_cigar_op(op)) {
        case BAM_CMATCH:
        case BAM_CEQUAL: {
            type = CIGAR_t::MATCH;
            break;
        }
        case BAM_CDIFF: {
            type = CIGAR_t::MISMATCH;
            break;
        }
        case BAM_CDEL: {
            type = CIGAR_t::DEL;
            break;
        }
        case BAM_CINS: {
            type = CIGAR_t::INS;
            break;
        }
        case BAM_CSOFT_CLIP: {
            type = CIGAR_t::SOFT;
            break;
        }
        case BAM_CHARD_CLIP: {
            type = CIGAR_t::HARD;
            break;
        }
        case BAM_CREF_SKIP: {
            type = CIGAR_t::SKIP;
            break;
        }
        case BAM_CPAD: {
            type = CIGAR_t::PAD;
            break;
        }
        default: {
            type = CIGAR_t::UNKNOWN;
            break;
        }
    }

    return CIGAR_op(type, length);

}


static inline std::span<const uint32_t> _get_cigar_str(bam1_t* _hts_aln) {

    uint32_t* raw_str = bam_get_cigar(_hts_aln);
    if (raw_str == nullptr) {
        throw std::runtime_error("Error retrieving the CIGAR string.");
    }
    hts_pos_t length = _hts_aln->core.n_cigar;

    return std::span<const uint32_t>(raw_str, raw_str + length);

}


static inline std::string _get_md_tag(bam1_t* _hts_aln) {

    uint8_t* md_ptr = bam_aux_get(_hts_aln, "MD");
    if (md_ptr == nullptr) {
        throw std::runtime_error("Failed to retrieve the MD tag.");
    }
    char* tag = bam_aux2Z(md_ptr);

    return std::string(tag, strlen(tag));

}


static inline std::string _cigar_hts_format(CIGAR_t type) {

    switch (type) {
        case CIGAR_t::MATCH: { return "M"; }
        case CIGAR_t::INS:   { return "I"; }
        case CIGAR_t::DEL:   { return "D"; }
        default:             { return "?"; }
    }

}





//
// Mapping and PHRED quality tools
//





template <typename dtype>
static inline std::vector<dtype> _get_mask(
    std::span<const qual_t> quality,
    qual_t min,
    int64_t window
) {

    int64_t length = quality.size();
    int64_t bad    = 0;

    window = std::min(window, length);
    int64_t start, stop;
    if (length < 2 * window + 1) {
        start = length - window;
        stop  = window + 1;
    } else {
        start = window + 1;
        stop  = length - window;
    }

    std::vector<dtype> mask(length + 1, 0);

    // 1. Count the bad bases in the starting window region
    for (int64_t ix = 0; ix < window; ix++) {
        bad += (quality[ix] < min);
    }

    // 2. Increment bad bases at the 3' end as the window region grows
    for (int64_t ix = 0; ix < start; ix++) {
        bad += (quality[ix + window] < min);
        mask[ix] = (bad == 0);
    }

    // 3. Increment and decrement bad bases as the window region shifts
    if (length < 2 * window + 1) {
        for (int64_t ix = start; ix < stop; ix++) {
            mask[ix] = (bad == 0);
        }
    } else {
        for (int64_t ix = start; ix < stop; ix++) {
            bad -= (quality[ix - window - 1] < min);
            bad += (quality[ix + window] < min);
            mask[ix] = (bad == 0);
        }
    }

    // 4. Decrement bad bases at the 5' end as the window region shrinks
    for (int64_t ix = stop; ix < length + 1; ix++) {
        bad -= (quality[ix - window - 1] < min);
        mask[ix] = (bad == 0);
    }

    return mask;

}


static inline std::string _phred_to_str(const std::vector<qual_t>& qualities) {

    int64_t _length = qualities.size();
    std::string str(_length, '\0');

    for (int64_t ix = 0; ix < _length; ix++) {
        str[ix] = static_cast<char>(qualities[ix] + PHRED_OFFSET);
    }

    return str;

}


PHRED::PHRED(const std::vector<qual_t>& qualities) : _qualities(qualities) {}


qual_t PHRED::operator[](int64_t ix) const {

    return _qualities[ix];

}

std::string PHRED::str() const {

    return _phred_to_str(_qualities);

}


bool PHRED::check(int64_t ix, qual_t min, int64_t window) const {

    int64_t _length = static_cast<int64_t>(_qualities.size());
    int64_t _zero   = static_cast<int64_t>(0);

    int64_t jx_min = std::max(ix - window, _zero);
    int64_t jx_max = std::min(ix + window + 1, _length);
    bool quality = true;

    for (int64_t jx = jx_min; jx < jx_max; jx++) {
        quality &= (_qualities[jx] >= min);
        quality &= (_qualities[jx] < MISSING_MAPQ);
    }

    return quality;

}





//
// CIGAR and MD tags
//





static inline bool _is_digit(const std::string& tag, int64_t pos) {

    return pos < static_cast<int64_t>(tag.size()) && std::isdigit(tag[pos]);

}


static inline bool _is_alpha(const std::string& tag, int64_t pos) {

    return pos < static_cast<int64_t>(tag.size()) && std::isalpha(tag[pos]);

}


static inline bool _is_deletion(const std::string& tag, int64_t pos) {

    return pos < static_cast<int64_t>(tag.size()) && tag[pos] == MD_DEL;

}


static inline bool _is_null(const std::string& tag, int64_t pos) {

    return pos < static_cast<int64_t>(tag.size()) && tag[pos] == MD_NULL;

}


static inline void _skip_null(const std::string& tag, int64_t& pos) {

    if (_is_null(tag, pos)) { pos++; }

}


static inline int64_t _count_matches(const std::string& tag, int64_t& pos) {

    int64_t _matches = 0;
    while (_is_digit(tag, pos)) {
        _matches = _matches * 10 + (tag[pos] - '0');
        pos++;
    }

    return _matches;

}


static inline int64_t _count_mismatches(const std::string& tag, int64_t& pos) {

    int64_t _mismatches = 0;
    while (_is_alpha(tag, pos)) {
        _mismatches++;
        pos++;
        _skip_null(tag, pos);
    }

    return _mismatches;

}


static inline int64_t _count_deletions(const std::string& tag, int64_t& pos) {

    int64_t _deletions = 0;
    pos++;
    while (_is_alpha(tag, pos)) {
        _deletions++;
        pos++;
    }

    _skip_null(tag, pos);
    return _deletions;

}


static inline CIGAR_t _md_type(const std::string& tag, int64_t pos) {

    if (_is_digit(tag, pos) && !_is_null(tag, pos)) {
        return CIGAR_t::MATCH;
    } else if (_is_alpha(tag, pos)) {
        return CIGAR_t::MISMATCH;
    } else if (_is_deletion(tag, pos)) {
        return CIGAR_t::DEL;
    } else {
        return CIGAR_t::UNKNOWN;
    }

}


MD_tag::MD_tag(const std::string& tag) : _tag(tag) { _skip_null(_tag, pos); }


std::string MD_tag::str() const { return _tag; }


void MD_tag::update() {

    CIGAR_t type = _md_type(_tag, pos);
    int64_t count = 0;
    switch (type) {
        case CIGAR_t::MATCH: {
            count = _count_matches(_tag, pos);
            break;
        }
        case CIGAR_t::MISMATCH: {
            count = _count_mismatches(_tag, pos);
            break;
        }
        case CIGAR_t::DEL: {
            count = _count_deletions(_tag, pos);
            break;
        }
        default: {
            throw std::runtime_error("Unknown MD tag operation \"" + std::to_string(_tag[pos]) + "\" at tag position " + std::to_string(pos) + ". The tag is " + _tag + ".");
        }
    }

    curr = CIGAR_op(type, count);

}


CIGAR_op MD_tag::advance() {

    // If the current op is empty, then move to the next op
    if (curr.empty()) { update(); }

    // Return the current op and set it to be empty
    CIGAR_op _tmp = curr;
    curr = CIGAR_op();
    return _tmp;

}


CIGAR_op MD_tag::advance(int64_t max) {

    // If the current op is empty, then move to the next op
    if (curr.empty()) { update(); }

    // Split the op based on the requested maximum size;
    // return the first half and store the second half
    std::pair<CIGAR_op, CIGAR_op> pair = curr.split(max);
    curr = pair.second;
    return pair.first;

}


CIGAR_op::CIGAR_op(CIGAR_t type) : _type(type), _length(1) {}


CIGAR_op::CIGAR_op(CIGAR_t type, int64_t length) : _type(type), _length(length) {}


CIGAR_t CIGAR_op::type() const { return _type; }


int64_t CIGAR_op::length() const { return _length; }


int64_t CIGAR_op::rlength() const {

    switch (_type) {

        case CIGAR_t::MATCH:
        case CIGAR_t::MISMATCH:
        case CIGAR_t::DEL:
        case CIGAR_t::SKIP: {
            return _length;
        }

        case CIGAR_t::INS:
        case CIGAR_t::SOFT:
        case CIGAR_t::HARD:
        case CIGAR_t::PAD:
        case CIGAR_t::BACK:
        case CIGAR_t::UNKNOWN: {
            return 0;
        }

    }

}


int64_t CIGAR_op::qlength() const {

    switch (_type) {

        case CIGAR_t::MATCH:
        case CIGAR_t::MISMATCH:
        case CIGAR_t::INS:
        case CIGAR_t::SOFT: {
            return _length;
        }

        case CIGAR_t::DEL:
        case CIGAR_t::SKIP:
        case CIGAR_t::HARD:
        case CIGAR_t::PAD:
        case CIGAR_t::BACK:
        case CIGAR_t::UNKNOWN: {
            return 0;
        }

    }

}


std::string CIGAR_op::str() const {

    return std::to_string(_length) + _cigar_hts_format(_type);

}


void CIGAR_op::extend(int64_t val) {

    _length += val;

}


int64_t CIGAR_op::pos() const {

    return _pos;

}


bool CIGAR_op::empty() const {

    return _length == 0;

}


bool CIGAR_op::advance(int64_t offset) {

    if (_pos < _length + offset) {
        _pos++;
        return true;
    }
    return false;

}


std::pair<CIGAR_op, CIGAR_op> CIGAR_op::split(int64_t n) const {

    int64_t l1 = std::min(_length, n);
    int64_t l2 = _length - l1;

    return std::pair<CIGAR_op, CIGAR_op>(
        CIGAR_op(_type, l1),
        CIGAR_op(_type, l2)
    );

}


CIGAR_op CIGAR_op::match() const {
    return CIGAR_op(CIGAR_t::MATCH, _length);
}


CIGAR::CIGAR(int64_t size) {

    _str.reserve(size);

}


int64_t CIGAR::rlength() const {

    int64_t _rlength = 0;
    for (const auto& op : _str) {
        _rlength += op.rlength();
    }

    return _rlength;

}


std::string CIGAR::str() const {

    std::string _out_str;
    CIGAR_op prev(CIGAR_t::MATCH, 0);

    for (const auto& op: _str) {

        if (op.type() == CIGAR_t::MATCH || op.type() == CIGAR_t::MISMATCH) {
            prev.extend(op.length());
            continue;
        } 

        if (!prev.empty()) {
            _out_str += prev.str();
            prev = CIGAR_op(CIGAR_t::MATCH, 0);
        }
        _out_str += op.str();

    }

    if (!prev.empty()) {
        _out_str += prev.str();
    }

    return _out_str;

}


void CIGAR::extend(CIGAR_op op) {

    if (size() > 0 && _str.back().type() == op.type()) {
        _str.back().extend(op.length());
    } else {
        _str.push_back(op);
    }

}


void CIGAR::append(CIGAR_op op) {

    _str.push_back(op);

}


int64_t CIGAR::size() const {

    return _str.size();

}


int64_t CIGAR::bases() const {

    int64_t bases = 0;
    for (const auto& op : _str) {
        bases += op.length();
    }

    return bases;

}


CIGAR_op CIGAR::operator[](int64_t ix) const {

    return _str[ix];

}


CIGAR_op CIGAR::back() const {

    if (!_str.empty()) { return _str.back(); }
    return CIGAR_op();

}


auto CIGAR::begin() -> decltype(_str.begin()) { return _str.begin(); }


auto CIGAR::begin() const -> decltype(_str.begin()) { return _str.begin(); }


auto CIGAR::end() -> decltype(_str.end()) { return _str.end(); }


auto CIGAR::end() const -> decltype(_str.end()) { return _str.end(); }


static inline CIGAR _get_cigar(
    std::span<const uint32_t> hts_cigar,
    const std::string& md_tag_str
) {

    MD_tag md_tag(md_tag_str);

    // Reserve the minimum amount of space required
    CIGAR cigar(hts_cigar.size());

    for (const auto& hts_op : hts_cigar) {

        CIGAR_op op = _cigar_from_hts(hts_op);
        // If the CIGAR string indicates a match, we must use the MD tag
        // in order to find mismatches
        if (op.type() == CIGAR_t::MATCH) {
            int64_t remaining = op.length();
            while (remaining > 0) {
                CIGAR_op _md_cig = md_tag.advance(remaining);
                cigar.append(_md_cig);
                remaining -= _md_cig.length();
            }
        }
        // Else, the HTS CIGAR string contains all the information we need.
        // The MD tag must be advanced if it is a deletion to keep the two in sync.
        else {
            cigar.append(op);
            if (op.type() == CIGAR_t::DEL) {
                (void)md_tag.advance();
            }
        }
    }

    return cigar;

}





















static inline bam1_t* _open_aln() {

    bam1_t* _hts_aln = bam_init1();
    if (_hts_aln == nullptr) {
        throw std::runtime_error("Failed to allocate memory for a sequence aligment.");
    }

    // Set the empty alignment to be unaligned
    _hts_aln->core.tid = -1;
    return _hts_aln;

}


static inline void _close_aln(bam1_t*& _hts_aln) {

    if (_hts_aln != nullptr) {
        bam_destroy1(_hts_aln);
        _hts_aln = nullptr;
    }

}





//
// Printing overloads
//





std::ostream& operator<<(std::ostream& os, FileType fileType) {
    switch (fileType) {
        case FileType::SAM:  {
            os << "SAM";
            break;
        }
        case FileType::BAM:  {
            os << "BAM";
            break;
        }
        case FileType::CRAM: {
            os << "CRAM";
            break;
        }
    }
    return os;
}


std::ostream& operator<<(std::ostream& os, CIGAR_t cigar) {

    switch (cigar) {
        case CIGAR_t::UNKNOWN:  {
            os << "UNKNOWN";
            break;
        }
        case CIGAR_t::MATCH:    {
            os << "MATCH";
            break;
        }
        case CIGAR_t::MISMATCH: {
            os << "MISMATCH";
            break;
        }
        case CIGAR_t::DEL:      {
            os << "DEL";
            break;
        }
        case CIGAR_t::INS:      {
            os << "INS";
            break;
        }
        case CIGAR_t::SOFT:     {
            os << "SOFT";
            break;
        }
        case CIGAR_t::HARD:     {
            os << "HARD";
            break;
        }
        case CIGAR_t::SKIP:     {
            os << "SKIP";
            break;
        }
        case CIGAR_t::PAD:      {
            os << "PAD";
            break;
        }
        case CIGAR_t::BACK:     {
            os << "BACK";
            break;
        }
    }

    return os;

}


std::ostream& operator<<(std::ostream& os, CIGAR cigar) {

    for (const auto& op : cigar) {
        os << std::left << std::setw(9) << op.type() << std::right << std::setw(4) << op.length() << "\n";
    }

    return os;

}





//
// Reading and writing binary sequence data
//





static inline base_t _to_int(char base) {

    switch (base) {
        case A:   { return IX_A;   }
        case C:   { return IX_C;   }
        case G:   { return IX_G;   }
        case T:   { return IX_T;   }
        case U:   { return IX_U;   }
        default:  { return IX_UNK; }
    }

}


static inline char _to_char(base_t base) {

    switch (base) {

        case IX_A: { return A; }
        case IX_C: { return C; }
        case IX_G: { return G; }
        case IX_T: { return T; }
        default:   { return N; }

    }

}


std::string _to_str(const seq_t& sequence) {

    std::string str(sequence.size(), '\0');
    for (hts_pos_t ix = 0; ix < sequence.size(); ix++) {
        str[ix] = _to_char(sequence[ix]);
    }

    return str;

}


static inline int64_t _bytes_from_seq(int64_t n) {

    return (n + BASES_PER_BYTE - 1) / BASES_PER_BYTE;

}


static inline base_t _binary_mask(base_t nuc, int64_t i) {

    return nuc << ((3 - (i % BASES_PER_BYTE)) * 2);

}


static inline seq_t _to_binary(const std::string& sequence) {

    int64_t len   = sequence.length();
    int64_t bytes = _bytes_from_seq(len);
    seq_t binary(bytes);

    for (int64_t i = 0; i < len; i++) {

        base_t nuc = _to_int(sequence[i]);
        if (nuc == IX_UNK) {
            throw std::runtime_error("Invalid nucleotide \"" + std::to_string(sequence[i]) + "\" encountered.");
        }

        binary[i / BASES_PER_BYTE] |= _binary_mask(nuc, i);

    }

    return binary;

}


static inline base_t _from_binary(base_t binary, int64_t j) {

    return (binary >> ((3 - j) * 2)) & BIN_T;

}


static inline int64_t _mod_round_up(int64_t val, int64_t div) {

    int64_t out = val % div;
    if (out == 0) { out = div; }
    return out;

}


static inline seq_t _from_binary(const seq_t& binary, int64_t length) {

    int64_t bytes = binary.size();
    seq_t sequence(length);

    int64_t ix, jx;
    for (ix = 0; ix < bytes - 1; ix++) {
        for (jx = 0; jx < BASES_PER_BYTE; jx++) {
            sequence[BASES_PER_BYTE * ix + jx] = _from_binary(binary[ix], jx);
        }
    }

    int64_t end = _mod_round_up(length, BASES_PER_BYTE);
    for (jx = 0; jx < end; jx++) {
        sequence[BASES_PER_BYTE * ix + jx] = _from_binary(binary[ix], jx);
    }

    return sequence;

}


static inline seq_t _read_binary_sequence(std::ifstream& file, int64_t length, int64_t offset) {

    int64_t bytes = _bytes_from_seq(length);
    seq_t binary(bytes);

    file.seekg(offset);
    file.read(reinterpret_cast<char*>(binary.data()), bytes * sizeof(base_t));

    return _from_binary(binary, length);

}


static inline void _write_binary_sequence(std::ofstream& file, const std::string& sequence) {

    seq_t binary = _to_binary(sequence);
    file.write(reinterpret_cast<const char*>(binary.data()), binary.size() * sizeof(base_t));

}































static inline int64_t _get_header_length(BGZF* _bgzf_file) {

    int32_t length;
    _read_bgzf(_bgzf_file, &length, sizeof(length));
    return static_cast<int64_t>(length);

}

static inline bool _is_sorted_rec(const std::string& str) {

    int64_t pos = str.find(HEADER_SORT_KEY) + HEADER_SORT_KEY.length();
    std::string value = str.substr(pos);

    if (value == BGZF_SORTED) {
        return true;
    }
    return false;

}

static inline bool _is_sorted(const std::string& name) {

    BGZF* _bgzf_file = _open_bgzf(name);

    char buffer[BGZF_BUFFER];
    int64_t ix = 0;
    int64_t length = _get_header_length(_bgzf_file);
    std::string line = _read_bgzf_line(_bgzf_file, buffer, ix, length);

    _close_bgzf(_bgzf_file);
    return _is_sorted_rec(line);

}

static inline int64_t _length_from_rec(const std::string& str) {

    int64_t pos = str.find(HEADER_LEN_KEY) + HEADER_LEN_KEY.length();
    std::string value = str.substr(pos);
    return static_cast<int64_t>(std::stoi(value));

}

static inline bool _is_seq_rec(const std::string& line) {

    return line.starts_with("@SQ");

}

static inline int64_t _get_reference_lengths(BGZF* _bgzf_file) {

    char buffer[BGZF_BUFFER];
    int64_t length = _get_header_length(_bgzf_file);
    int64_t ix = 0;

    std::string line = _read_bgzf_line(_bgzf_file, buffer, ix, length);

    // Read header line by line
    while (ix < length) {
        int64_t _to_read = std::min(BGZF_BUFFER, length - ix);
        int64_t bytes = bgzf_read(_bgzf_file, buffer, _to_read);
        ix += bytes;
    }

    uint32_t references, name_len, ref_len;
    _read_bgzf(_bgzf_file, &references, sizeof(uint32_t));
    for (int32_t ix = 0; ix < references; ix++) {
        _read_bgzf(_bgzf_file, &name_len, sizeof(uint32_t));
        _read_bgzf(_bgzf_file, buffer, name_len);
        _read_bgzf(_bgzf_file, &ref_len, sizeof(uint32_t));
    }

    return references;

}


struct HeaderData {

    bool    sorted     = false;
    int64_t references = 0;
    int64_t ptr        = 0;

};


static inline HeaderData _read_header(BGZF* _bgzf_file, FileType type) {

    char buffer[BGZF_BUFFER];
    int64_t length = _get_header_length(_bgzf_file);
    int64_t ix = 0;
    HeaderData data;

    // The first line tells us whether the file is sorted
    std::string line = _read_bgzf_line(_bgzf_file, buffer, ix, length);
    data.sorted = _is_sorted_rec(line);

    // Read header line by line
    while (ix < length) {

        line = _read_bgzf_line(_bgzf_file, buffer, ix, length);
        if (_is_seq_rec(line)) { data.references++; }

    }

    if (type == FileType::BAM) {

        uint32_t references, name_len, ref_len;
        _read_bgzf(_bgzf_file, &references, sizeof(uint32_t));
        for (int32_t ix = 0; ix < references; ix++) {
            _read_bgzf(_bgzf_file, &name_len, sizeof(uint32_t));
            _read_bgzf(_bgzf_file, buffer, name_len);
            _read_bgzf(_bgzf_file, &ref_len, sizeof(uint32_t));
        }

    }

    data.ptr = bgzf_tell(_bgzf_file);
    return data;

}





//
// FASTA
//





static inline HeaderBlock _read_tiny_header_block(std::ifstream& file) {

    HeaderBlock block;
    file.read(reinterpret_cast<char*>(&block.length), sizeof(int64_t));
    file.read(reinterpret_cast<char*>(&block.sequences), sizeof(int64_t));
    return block;

}


static inline std::vector<HeaderBlock> _read_tiny_header(std::ifstream& file) {

    HeaderBlock block;
    std::vector<HeaderBlock> blocks;

    while (block.length != EOH) {

        block = _read_tiny_header_block(file);
        if (!block.empty()) {
            blocks.push_back(block);
        }

    }

    return blocks;

}


static inline void _write_tiny_header_block(std::ofstream& file, const HeaderBlock& block) {

    file.write(reinterpret_cast<const char*>(&block.length), sizeof(int64_t));
    file.write(reinterpret_cast<const char*>(&block.sequences), sizeof(int64_t));

}


static inline void _write_tiny_eoh(std::ofstream& file) {

    file.write(reinterpret_cast<const char*>(&EOH), sizeof(int64_t));

}


static inline void _write_tiny_header(std::ofstream& file, const std::vector<HeaderBlock>& blocks) {

    for (const auto& block : blocks) {
        if (!block.empty()) {
            _write_tiny_header_block(file, block);
        }
    }
    _write_tiny_eoh(file);

}


static inline std::vector<HeaderBlock> _fasta_to_tiny_header_blocks(const std::string& name) {

    OldSchool::FASTA fasta(name);
    HeaderBlock block;
    std::vector<HeaderBlock> blocks;

    std::string sequence = fasta.next();
    while (!sequence.empty()) {

        int64_t length = sequence.length();
        if (length == block.length) {
            block.sequences++;
        } else {
            blocks.push_back(block);
            block.length = length;
            block.sequences = 1;
        }
        sequence = fasta.next();

    }

    if (!block.empty()) {
        blocks.push_back(block);
    }

    return blocks;

}


static inline void _fasta_to_tiny(const std::string& fasta_name, const std::string& tiny_name) {

    if (fasta_name == tiny_name) {
        throw std::runtime_error("The input and output files must be different.");
    }
    _throw_if_not_exists(fasta_name);
    _throw_if_exists(tiny_name);

    std::ofstream tiny(tiny_name);
    std::vector<HeaderBlock> blocks = _fasta_to_tiny_header_blocks(fasta_name);
    _write_tiny_header(tiny, blocks);

    OldSchool::FASTA fasta(fasta_name);
    std::string sequence = fasta.next();
    if (sequence.empty()) { return; }

    _write_binary_sequence(tiny, sequence);
    while (!sequence.empty()) {
        sequence = fasta.next();
        _write_binary_sequence(tiny, sequence);
    }

    tiny.close();

}





bool HeaderBlock::empty() const {

    return (length <= 0 || sequences == 0);

}


FASTA::FASTA(const std::string& name) {

    _throw_if_not_exists(name);

    _fasta_name = name;
    _name       = name + ".tiny";

    if (!_exists(_name)) {
        _fasta_to_tiny(_fasta_name, _name);
    }

    _file   = std::ifstream(_name);
    _blocks = _read_tiny_header(_file);
    _offset = (2 * _blocks.size() + 1) * sizeof(int64_t);

}


FASTA::FASTA(FASTA&& other) noexcept
    : _name(std::move(other._name)),
      _file(std::move(other._file)),
      _blocks(std::move(other._blocks)),
      _offset(other._offset) {}


FASTA& FASTA::operator=(FASTA&& other) noexcept {

    if (this != &other) {
        _name = other._name;
        _file = std::move(other._file);
        _blocks = other._blocks;
        _offset = other._offset;
    }

    return *this;

}


FASTA::FASTA(const FASTA& other)
    : _name(other._name),
      _file(other._name),
      _blocks(other._blocks),
      _offset(other._offset) {}


std::string FASTA::name() const {

    return _name;

}


int64_t FASTA::size() const {

    int64_t _sequences = 0;
    for (const auto& block : _blocks) {
        _sequences += block.sequences;
    }

    return _sequences;

}


int64_t FASTA::length(int64_t ix) const {

    int64_t count = 0;
    for (const auto& block : _blocks) {

        if (count + block.sequences > ix) {
            return block.length;
        }
        count += block.sequences;

    }

    return -1;

}


seq_t FASTA::sequence(int64_t ix) {

    int64_t tmp = 0, count = 0, length = -1;
    for (const auto& block : _blocks) {

        int64_t bytes = _bytes_from_seq(block.length);

        if (count + block.sequences > ix) {
            tmp += static_cast<int64_t>(ix - count) * bytes;
            length = block.length;
            break;
        } else {
            tmp += block.sequences * bytes;
            count += block.sequences;
        }

    }

    // Account for the header
    int64_t offset = tmp * sizeof(base_t) + _offset;
    return _read_binary_sequence(_file, length, offset);

}


int64_t FASTA::longest() const {

    int64_t _longest = 0;
    for (const auto& block : _blocks) {

        if (block.length > _longest) {
            _longest = block.length;
        }

    }

    return _longest;

}









static inline void _build_tiny_index(
    BGZF* _bgzf_file,
    const std::string& filename,
    int32_t references
) {

    _throw_if_exists(filename);
    std::ofstream outfile(filename);

    bam1_t* _hts_aln = _open_aln();
    IndexBlock block;
    int64_t unaligned = 0;

    // The beginning of the non-header section
    block.tell(_bgzf_file);
    outfile.seekp(0);
    block.write_ptr(outfile);

    int32_t tid = 0, curr;
    while (bam_read1(_bgzf_file, _hts_aln) > 0) {

        curr = _hts_aln->core.tid;
        if (curr > tid) {

            block.write_reads(outfile);
            block.reads = 0;

            // Write empty blocks for all references which did not appear in the file
            for (int32_t jx = tid + 1; jx < curr; jx++) {
                block.write_ptr(outfile);
                block.write_reads(outfile);
            }

            tid = curr;
            block.write_ptr(outfile);

        } else if (curr == -1) {

            unaligned++;

        }

        block.reads++;
        block.tell(_bgzf_file);

    }

    block.write_reads(outfile);
    block.reads = 0;

    for (int32_t jx = tid + 1; jx < references; jx++) {
        block.write_ptr(outfile);
        block.write_reads(outfile);
    }

    outfile.write(reinterpret_cast<char*>(&unaligned), sizeof(int64_t));
    outfile.close();

}




//
// Index
//





bool IndexBlock::empty() const {

    return (ptr == -1 || reads == 0);

}

void IndexBlock::tell(BGZF* _hts_bgzf) {

    ptr = bgzf_tell(_hts_bgzf);

}

void IndexBlock::write_ptr(std::ofstream& file) {

    file.write(reinterpret_cast<char*>(&ptr), sizeof(int64_t));

}

void IndexBlock::write_reads(std::ofstream& file) {

    file.write(reinterpret_cast<char*>(&reads), sizeof(int64_t));

}


Index::Index(const std::string& filename, int64_t references) : _name(filename), _file(filename), _references(references) {}


IndexBlock Index::read(int64_t ix) {

    if (ix >= _references) {
        throw std::runtime_error("The index " + std::to_string(ix) + " is too large for the file \"" + _name + "\".");
    }

    IndexBlock block;

    _file.seekg(ix * sizeof(IndexBlock));
    _file.read(reinterpret_cast<char*>(&block), sizeof(IndexBlock));

    return block;

}


int64_t Index::aligned() {

    int64_t reads = 0;
    for (int64_t ix = 0; ix < _references; ix++) {
        reads += read(ix).reads;
    }

    return reads;

}


int64_t Index::unaligned() {

    int64_t unaligned;

    _file.seekg(_references * sizeof(IndexBlock));
    _file.read(reinterpret_cast<char*>(&unaligned), sizeof(int64_t));

    return unaligned;

}





//
// Alignment
//






Alignment::Alignment(BGZF* _hts_bgzf, bam1_t* _hts_aln, int64_t reads) : _hts_bgzf(_hts_bgzf), _hts_aln(_hts_aln), _reads(reads) {}


Alignment::Alignment(Alignment&& other) noexcept
    : _hts_bgzf(other._hts_bgzf),
      _hts_aln(other._hts_aln),
      _reads(other._reads),
      _curr(other._curr) {

    other._hts_bgzf = nullptr;
    other._hts_aln  = nullptr;

}


Alignment& Alignment::operator=(Alignment&& other) noexcept {

    if (this != &other) {
        _hts_bgzf = other._hts_bgzf;
        _hts_aln  = other._hts_aln;
        _reads = other._reads;
        _curr = other._curr;
        other._hts_bgzf = nullptr;
        other._hts_aln  = nullptr;
    }

    return *this;

}


bool Alignment::next() {

    if (_curr < _reads) {

        if (bam_read1(_hts_bgzf, _hts_aln) < 0) {
            throw std::runtime_error("The alignment ended prematurely.");
        }
        _curr++;
        return true;

    }

    return false;

}


bool Alignment::aligned() const {

    return _hts_aln->core.tid >= 0;

}


base_t Alignment::base(int64_t ix) const {

    uint8_t* _query = bam_get_seq(_hts_aln);
    uint8_t _hts_base = bam_seqi(_query, ix);

    return _base_from_hts(_hts_base);

}


int64_t Alignment::length() const {

    return _hts_aln->core.l_qseq;

}


int64_t Alignment::offset() const {

    return _hts_aln->core.pos;

}


qual_t Alignment::mapq() const {

    return _hts_aln->core.qual;

}


std::span<const qual_t> Alignment::phred() const {

    qual_t* _quality = bam_get_qual(_hts_aln);
    return std::span<const qual_t>(_quality, _quality + length());

}


template <typename dtype>
std::vector<dtype> Alignment::mask(qual_t min, int64_t window) const {

    return _get_mask<dtype>(phred(), min, window);

}


CIGAR Alignment::cigar() const {

    std::string _md_tag = _get_md_tag(_hts_aln);
    std::span<const uint32_t> _cigar = _get_cigar_str(_hts_aln);

    return _get_cigar(_cigar, _md_tag);

}





//
// File
//





File::File(const std::string& name) : _name(name) {

    _hts_bgzf = _open_bgzf(name);
    _hts_aln  = _open_aln();
    _type     = _get_filetype(_hts_bgzf);

    HeaderData _data = _read_header(_hts_bgzf, _type);
    _sorted     = _data.sorted;
    _references = _data.references;
    _body_ptr   = _data.ptr;

    std::string _index_name = _name + ".tinyindex";
    if (!std::filesystem::exists(_index_name)) {
        _build_tiny_index(_hts_bgzf, _index_name, _references);
    }
    _index = Index(_index_name, _references);

}


File::File(File&& other) noexcept
    : 
      _name(std::move(other._name)),
      _hts_bgzf(other._hts_bgzf),
      _hts_aln(other._hts_aln),
      _type(other._type),
      _sorted(other._sorted),
      _references(other._references),
      _body_ptr(other._body_ptr),
      _index(std::move(other._index))
{

    other._hts_bgzf   = nullptr;
    other._hts_aln  = nullptr;

}


File::~File() { 

    _close_aln(_hts_aln);
    _close_bgzf(_hts_bgzf);

}


std::string File::name() const { return _name; }


BGZF* File::ptr() const { return _hts_bgzf; }


FileType File::type() const { return _type; }


Alignment File::alignment(int64_t ix, bool seek) {

    IndexBlock block = _index.read(ix);
    if (seek) { _seek_bgzf(_hts_bgzf, block.ptr); }
    return Alignment(_hts_bgzf, _hts_aln, block.reads);

}


int64_t File::size() const { return _references; }


int64_t File::aligned()   { return _index.aligned(); }


int64_t File::unaligned() { return _index.unaligned(); }





//
// FileGroup
//





FileGroup::FileGroup(const std::vector<std::string>& filenames) {

    _throw_if_has_duplicate_paths(filenames);
    for (auto& name : filenames) { _group.emplace_back(name); }

}


int64_t FileGroup::size() const { return _group.size(); }


int64_t FileGroup::aligned() {

    int64_t _reads = 0;
    for (auto& file : _group) {
        _reads += file.aligned();
    }

    return _reads;

}


int64_t FileGroup::unaligned() {

    int64_t _reads = 0;
    for (auto& file : _group) {
        _reads += file.unaligned();
    }

    return _reads;

}


int64_t FileGroup::references() const {

    if (size() > 0) {
        return _group[0].size();
    } else {
        return 0;
    }

}


auto FileGroup::begin() -> decltype(_group.begin()) { return _group.begin(); }


auto FileGroup::begin() const -> decltype(_group.begin()) { return _group.begin(); }


auto FileGroup::end() -> decltype(_group.end()) { return _group.end(); }


auto FileGroup::end() const -> decltype(_group.end()) { return _group.end(); }





//
// Template specialization
//





template std::vector<float> Alignment::mask<float>(qual_t, int64_t) const;





} // namespace TinyHTS
