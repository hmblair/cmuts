#include "hts.hpp"
#include "mpi.hpp"
#include <random>
#include <stdexcept>

namespace HTS {

void __disable_logging() {

    hts_set_log_level(htsLogLevel::HTS_LOG_OFF);

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

static inline FileType _filetype_from_hts(int fmt) {

    switch (fmt) {

        case sam:  { return FileType::SAM;  }
        case bam:  { return FileType::BAM;  }
        case cram: { return FileType::CRAM; }
        default:   {
            throw std::runtime_error("Invalid HTS format.");
        }

    }

}

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

static inline bool _is_mod(CIGAR_t type) {

    switch (type) {

        case CIGAR_t::MISMATCH:
        case CIGAR_t::INS:
        case CIGAR_t::DEL:     {
            return true;
        }

        case CIGAR_t::MATCH:
        case CIGAR_t::SOFT:
        case CIGAR_t::HARD:
        case CIGAR_t::SKIP:
        case CIGAR_t::PAD:
        case CIGAR_t::BACK:
        case CIGAR_t::UNKNOWN: {
            return false;
        }

    }

}

static inline bool _is_digit(const std::string& tag, hts_pos_t pos) {
    return pos < tag.size() && std::isdigit(tag[pos]);
}

static inline bool _is_alpha(const std::string& tag, hts_pos_t pos) {
    return pos < tag.size() && std::isalpha(tag[pos]);
}

static inline bool _is_deletion(const std::string& tag, hts_pos_t pos) {
    return pos < tag.size() && tag[pos] == MD_DEL;
}

static inline bool _is_null(const std::string& tag, hts_pos_t pos) {
    return pos < tag.size() && tag[pos] == MD_NULL;
}

static inline void _skip_null(const std::string& tag, hts_pos_t& pos) {
    if (_is_null(tag, pos)) { pos++; }
}

static inline hts_pos_t _count_matches(const std::string& tag, hts_pos_t& pos) {
    hts_pos_t _matches = 0;
    while (_is_digit(tag, pos)) {
        _matches = _matches * 10 + (tag[pos] - '0');
        pos++;
    }
    return _matches;
}

static inline hts_pos_t _count_mismatches(const std::string& tag, hts_pos_t& pos) {
    hts_pos_t _mismatches = 0;
    while (_is_alpha(tag, pos)) {
        _mismatches++;
        pos++;
        _skip_null(tag, pos);
    }
    return _mismatches;
}

static inline hts_pos_t _count_deletions(const std::string& tag, hts_pos_t& pos) {
    hts_pos_t _deletions = 0;
    pos++;
    while (_is_alpha(tag, pos)) {
        _deletions++;
        pos++;
    }
    _skip_null(tag, pos);
    return _deletions;
}

static inline CIGAR_t _md_type(const std::string& tag, hts_pos_t pos) {
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


MD::MD(const std::string& tag) : _tag(tag) {
    _skip_null(_tag, pos);
}

std::string MD::tag() const {
    return _tag;
}

void MD::update() {

    CIGAR_t type = _md_type(_tag, pos);
    hts_pos_t count = 0;
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

CIGAR_op MD::advance(hts_pos_t max) {

    // If the current op is empty, then move to the next op
    if (curr.empty()) { update(); }

    // Split the op based on the requested maximum size;
    // return the first half and store the second half
    std::pair<CIGAR_op, CIGAR_op> pair = curr.split(max);
    curr = pair.second;
    return pair.first;

}


// Helper functions

static inline void _throw_if_exists(const std::string& filename) {
    if (std::filesystem::exists(filename)) {
        throw std::runtime_error("The file \"" + filename + "\" already exists.");
    }
}

static inline void _throw_if_not_exists(const std::string& filename) {
    if (!std::filesystem::exists(filename)) {
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
    if (parent.empty()) {
        return stem;
    }
    return parent + "/" + stem;
}

static inline bool _has_duplicates(const std::vector<std::string>& paths) {

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

static inline void _throw_if_has_duplicates(const std::vector<std::string>& paths) {
    if (_has_duplicates(paths)) {
        throw std::runtime_error("Duplicate input paths (modulo extensions) are not allowed.");
    }
}

static inline base_t _as_int(char base) {
    switch (base) {
        case A:   { return IX_A;   }
        case C:   { return IX_C;   }
        case G:   { return IX_G;   }
        case T:   { return IX_T;   }
        case U:   { return IX_U;   }
        default:  { return IX_UNK; }
    }
}

static inline seq_t _as_int(
    const std::string& sequence
) {
    int64_t length = sequence.length();
    seq_t result(length);
    for (int64_t i = 0; i < length; i++) {
        result[i] = _as_int(sequence[i]);
    }
    return result;
}

static inline void _as_int(
    const std::string& sequence,
    view_t<base_t, DS_DIMS> buffer
) {
    for (int64_t i = 0; i < sequence.length(); i++) {
        buffer(i) = _as_int(sequence[i]);
    }
}

static inline faidx_t* _open_fai(const std::string& name) {
    faidx_t* _hts_fai = fai_load(name.c_str());
    if (_hts_fai == nullptr) {
        throw std::runtime_error("Failed to load the .fai index associated with \"" + name + "\".");
    }
    return _hts_fai;
}

static inline void _close_fai(faidx_t*& _hts_fai) {
    if (_hts_fai != nullptr) {
        fai_destroy(_hts_fai);
        _hts_fai = nullptr;
    }
}

static inline const char* _get_name(const FASTA& fasta, hts_pos_t ix) {
    const char *name = faidx_iseq(fasta.ptr(), ix);
    if (name == nullptr) {
        throw std::runtime_error(
            "There was an error retrieving the name of sequence " + std::to_string(ix) + " in \"" + fasta.name() + "\"."
        );
    }
    return name;
}

static inline hts_pos_t _length_from_name(const FASTA& fasta, const char* name) {
    int len = 0;
    char *sequence = fai_fetch(fasta.ptr(), name, &len);
    if (sequence == nullptr || len < 0) {
        throw std::runtime_error(
            "There was an error fetching a sequence from \"" + fasta.name() + "\"."
        );
    }
    return len;
}

static inline hts_pos_t _get_length(const FASTA& fasta, hts_pos_t ix) {
    const char* name = _get_name(fasta, ix);
    return _length_from_name(fasta, name);
}

static inline std::string _sequence_from_name(const FASTA& fasta, const char* name) {
    int len = 0;
    char *sequence = fai_fetch(fasta.ptr(), name, &len);
    if (sequence == nullptr || len < 0) {
        throw std::runtime_error(
            "There was an error fetching a sequence from \"" + fasta.name() + "\"."
        );
    }
    std::string result(sequence, len);
    free(sequence);
    return result;
}

static inline std::string _get_sequence(const FASTA& fasta, hts_pos_t ix) {
    const char* name = _get_name(fasta, ix);
    return _sequence_from_name(fasta, name);
}

static inline void _verify_index(const FASTA& fasta, hts_pos_t ix) {
    if (ix < 0 || ix >= fasta.size()) {
        throw std::runtime_error(
            "The index (" + std::to_string(ix) + 
            ") is not valid given the number of sequences in the FASTA file \"" + 
            fasta.name() + "\" (" + std::to_string(fasta.size()) + ")."
        );
    }
}


static inline bool _has_index(FileType type, const std::string& name) {
    switch (type) {
        case FileType::SAM:
        case FileType::BAM: {
            return std::filesystem::is_regular_file(name + BAM_INDEX);
        }
        case FileType::CRAM: {
            return std::filesystem::is_regular_file(name + CRAM_INDEX);
        }
    }
}

static inline void _build_index(FileType type, const std::string& name) {

    if (_has_index(type, name)) {
        return;
    }

    if (sam_index_build(name.c_str(), 0) < 0) {
        throw std::runtime_error("Failed to build the index for the file \"" + name + "\".");
    }

}

static inline htsFile* _open_file(const std::string& name) {

    htsFile* _hts_file = hts_open(name.c_str(), "r");
    if (_hts_file == nullptr) {
        throw std::runtime_error("Failed to open the file \"" + name + "\".");
    }

    return _hts_file;

}

static inline void _close_file(htsFile*& _hts_file) {
    if (_hts_file != nullptr) {
        sam_close(_hts_file);
        _hts_file = nullptr;
    }
}


static inline bam_hdr_t* _open_header(htsFile* _hts_file, const std::string& name) {

    bam_hdr_t* _hts_header = sam_hdr_read(_hts_file);
    if (_hts_header == nullptr) {
        throw std::runtime_error("Failed to read the header for \"" + name + "\".");
    }

    return _hts_header;

}

static inline void _close_header(bam_hdr_t*& _hts_header) {

    if (_hts_header != nullptr) {
        bam_hdr_destroy(_hts_header);
        _hts_header = nullptr;
    }

}

static inline bool _is_sorted(bam_hdr_t* _hts_header) {

    if (_hts_header == nullptr) {
        throw std::runtime_error("The header must be opened to check if the file is sorted.");
    }

    kstring_t ks = {0, 0, NULL};
    return sam_hdr_find_line_id(_hts_header, "HD", "SO", "coordinate", &ks) == 0;

}

static inline void _throw_if_not_sorted(bam_hdr_t* _hts_header, const std::string& name) {

    if (!_is_sorted(_hts_header)) {
        throw std::runtime_error("The input file \"" + name + "\" is not sorted.");
    }

}

static inline void _sort_core(const std::string& name, bam_hdr_t* _hts_header, const MPI::Manager& mpi) {

    if (mpi.root() && !_is_sorted(_hts_header)) {

        std::string unsorted_name = _path(name) + "_UNSORTED.bam";
        _safe_move(name, unsorted_name);

        std::string command = "samtools sort -o " + name + " " + unsorted_name;
        if (system(command.c_str()) != 0) {
            _safe_move(unsorted_name, name);
            throw std::runtime_error("There was an error sorting the file \"" + name + "\".");
        }

    }

    mpi.barrier();

}

static inline void _sort(const std::string& name, const MPI::Manager& mpi) {

    htsFile* _hts_file = _open_file(name);
    bam_hdr_t* _hts_header = _open_header(_hts_file, name);
    _sort_core(name, _hts_header, mpi);
    _close_file(_hts_file);
    _close_header(_hts_header);

}

static inline hts_idx_t* _open_index(const File& file) {

    _throw_if_not_sorted(file.header(), file.name());

    if (file.mpi().root()) {
        _build_index(file.type(), file.name());
    }
    file.mpi().barrier();

    hts_idx_t* _hts_index = sam_index_load(file.ptr(), file.name().c_str());
    if (_hts_index == nullptr) {
        throw std::runtime_error("Failed to open the index associated with \"" + file.name() + "\".");
    }

    return _hts_index;

}

static inline void _close_index(hts_idx_t*& _hts_index) {
    if (_hts_index != nullptr) {
        hts_idx_destroy(_hts_index);
        _hts_index = nullptr;
    }
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

static inline int64_t _sam_bam_aligned(hts_idx_t* _hts_index, hts_pos_t ix) {

    uint64_t mapped = 0, unmapped = 0;
    bool success = hts_idx_get_stat(_hts_index, ix, &mapped, &unmapped) >= 0;

    if (success) {
        return mapped;
    }
    return 0;

}

static inline int64_t _sam_bam_aligned(hts_idx_t* _hts_index) {

    int64_t aligned = 0;
    hts_pos_t size = hts_idx_nseq(_hts_index);

    for (int64_t ix = 0; ix < size; ix++) {
        aligned += _sam_bam_aligned(_hts_index, ix);
    }

    return aligned;

}

static inline int64_t _cram_aligned(htsFile* _hts_file, bam_hdr_t* _hts_header) {

    bam1_t* _hts_aln = _open_aln();
    int64_t _reads = 0;

    while (sam_read1(_hts_file, _hts_header, _hts_aln) > 0) {
        if (_hts_aln->core.tid >= 0) {
            _reads++;
        }
    }

    return _reads;

}

static inline int64_t _unaligned(hts_idx_t* _hts_index) {

    return static_cast<int64_t>(
        hts_idx_get_n_no_coor(_hts_index)
    );

}

static inline void _set_reference(const File& file, const std::string& fasta) {

    if (hts_set_fai_filename(file.ptr(), fasta.c_str()) < 0) {
        throw std::runtime_error("Failed to attach the reference FASTA to the file \"" + file.name() + "\".");
    }

}

static inline hts_itr_t* _open_iter(const File& file, hts_pos_t ix) {

    hts_itr_t* _hts_iter = sam_itr_queryi(file.index(), ix, 0, HTS_POS_MAX);
    if (_hts_iter == nullptr) {
        throw std::runtime_error("Failed to get the iterator for sequence " + std::to_string(ix) + " in \"" + file.name() + "\".");
    }

    return _hts_iter;

}

static inline void _close_iter(hts_itr_t*& _hts_iter) {

    if (_hts_iter != nullptr) {
        hts_itr_destroy(_hts_iter);
        _hts_iter = nullptr;
    }

}

void _throw_if_bad_lengths(const File& file) {
    const FASTA& fasta = file.fasta();
    if (file.size() != fasta.size()) {
        throw std::runtime_error("The number of references in the header for the file \"" + file.name() + "\" (" + std::to_string(file.size()) + ") does not match the number of sequences in the FASTA (" + std::to_string(fasta.size()) + ").");
    }
    for (hts_pos_t ix = 0; ix < file.size(); ix++) {
        if (file.length(ix) != fasta.length(ix)) {
            throw std::runtime_error("The sequence lengths in the header for the file \"" + file.name() + "\" and the FASTA file do not match.");
        }
    }
}

// Convert from HTS format (1,2,4,8) to standard
// format (0,1,2,3)
static inline base_t _base_from_hts(uint8_t base) {

    switch (base) {
        case HTS_A: { return IX_A;   }
        case HTS_C: { return IX_C;   }
        case HTS_G: { return IX_G;   }
        case HTS_T: { return IX_T;   }
        default:    { return IX_UNK; }
    }

}

static inline std::span<const uint32_t> _get_cigar_str(bam1_t* _hts_aln) {

    uint32_t* raw_str = bam_get_cigar(_hts_aln);
    if (raw_str == nullptr) {
        throw std::runtime_error("Error retrieve the CIGAR string.");
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

template <typename dtype>
static inline std::vector<dtype> _get_mask(
    std::span<const uint8_t> quality,
    uint8_t min,
    hts_pos_t window
) {

    hts_pos_t length = quality.size();
    hts_pos_t bad    = 0;

    window = std::min(window, length);
    hts_pos_t start, stop;
    if (length < 2 * window + 1) {
        start = length - window;
        stop  = window + 1;
    } else {
        start = window + 1;
        stop  = length - window;
    }

    std::vector<dtype> mask(length + 1, 0);

    // 1. Count the bad bases in the starting window region
    for (hts_pos_t ix = 0; ix < window; ix++) {
        bad += (quality[ix] < min);
    }

    // 2. Increment bad bases at the 3' end as the window region grows
    for (hts_pos_t ix = 0; ix < start; ix++) {
        bad += (quality[ix + window] < min);
        mask[ix] = (bad == 0);
    }

    // 3. Increment and decrement bad bases as the window region shifts
    if (length < 2 * window + 1) {
        for (hts_pos_t ix = start; ix < stop; ix++) {
            mask[ix] = (bad == 0);
        }
    } else {
        for (hts_pos_t ix = start; ix < stop; ix++) {
            bad -= (quality[ix - window - 1] < min);
            bad += (quality[ix + window] < min);
            mask[ix] = (bad == 0);
        }
    }

    // 4. Decrement bad bases at the 5' end as the window region shrinks
    for (hts_pos_t ix = stop; ix < length + 1; ix++) {
        bad -= (quality[ix - window - 1] < min);
        mask[ix] = (bad == 0);
    }

    return mask;

}

static inline CIGAR_op _cigar_from_hts(uint32_t op) {

    hts_pos_t length = bam_cigar_oplen(op);
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

static inline CIGAR _get_cigar(
    std::span<const uint32_t> hts_cigar,
    const std::string& md_tag
) {

    MD md(md_tag);

    // Reserve the minimum amount of space required
    CIGAR cigar(hts_cigar.size());

    for (const auto& hts_op : hts_cigar) {

        CIGAR_op op = _cigar_from_hts(hts_op);
        // If the CIGAR string indicates a match, we must use the MD tag
        // in order to find mismatches
        if (op.type() == CIGAR_t::MATCH) {
            hts_pos_t remaining = op.length();
            while (remaining > 0) {
                CIGAR_op _md_cig = md.advance(remaining);
                cigar.append(_md_cig);
                remaining -= _md_cig.length();
            }
        }
        // Else, the HTS CIGAR string contains all the information we need.
        // The MD tag must be advanced if it is a deletion to keep the two in sync.
        else {
            cigar.append(op);
            if (op.type() == CIGAR_t::DEL) {
                (void)md.advance();
            }
        }
    }

    return cigar;

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

static inline std::string _default_fasta_name(size_t n) {

    return ">ref" + std::to_string(n);

}

void _default_write_to_fasta(
    const std::string& filename,
    const std::vector<std::string>& sequences
) {

    _throw_if_exists(filename);

    std::ofstream file(filename);
    for (size_t ix = 0; ix < sequences.size(); ix++) {
        file << _default_fasta_name(ix) << "\n";
        file << sequences[ix] << "\n";
    }

}

static inline std::string _phred_as_str(const std::vector<uint8_t>& qualities) {

    size_t _length = qualities.size();
    std::string str(_length, '\0');

    for (int64_t ix = 0; ix < _length; ix++) {
        str[ix] = static_cast<char>(qualities[ix] + PHRED_OFFSET);
    }

    return str;

}


//
// PHRED
//



PHRED::PHRED(const std::vector<uint8_t>& qualities) : _qualities(qualities) {}

uint8_t PHRED::operator[](size_t ix) const {

    return _qualities[ix];

}

std::string PHRED::str() const {

    return _phred_as_str(_qualities);

}

bool PHRED::check(size_t ix, uint8_t min, size_t window) const {

    int64_t _ix     = static_cast<int64_t>(ix);
    int64_t _window = static_cast<int64_t>(window);
    int64_t _length = static_cast<int64_t>(_qualities.size());
    int64_t _zero   = static_cast<int64_t>(0);

    int64_t jx_min = std::max(_ix - _window, _zero);
    int64_t jx_max = std::min(_ix + _window + 1, _length);
    bool quality = true;

    for (int64_t jx = jx_min; jx < jx_max; jx++) {
        quality &= (_qualities[jx] >= min);
        quality &= (_qualities[jx] < MISSING_MAPQ);
    }

    return quality;

}



//
// FASTA
//



template <typename dtype>
static inline HDF5::Memspace<dtype, DS_DIMS> __memspace(
    const HTS::FASTA& fasta,
    HDF5::File& hdf5
) {
    std::vector<size_t> dims = {fasta.size(), static_cast<size_t>(fasta.max_length())};
    return hdf5.memspace<dtype, DS_DIMS>(dims, SEQUENCE_DS);
}

template <typename dtype>
static inline void _fasta_to_hdf5(
    const HTS::FASTA& fasta,
    HDF5::File& hdf5,
    const MPI::Manager& mpi
) {

    HDF5::Memspace<dtype, DS_DIMS> memspace = __memspace<dtype>(fasta, hdf5);
    MPI::Chunk chunk = mpi.chunk(
        hdf5.chunk_size(),
        fasta.size()
    );

    for (int64_t ix = chunk.low; ix < chunk.high; ix+= chunk.step) {
        for (int64_t jx = 0; jx < chunk.size; jx++) {
            if (ix + jx < fasta.size()) {
                view_t<dtype, DS_DIMS> arr = memspace.view(jx);
                fasta.as_int(ix + jx, arr);
            }
        }
        memspace.safe_write(ix);
        memspace.clear();
    }

}

FASTA::FASTA(const std::string& filename) : _name(filename) {

    _throw_if_not_exists(filename);
    _hts_fai = _open_fai(filename);
    _size    = faidx_nseq(_hts_fai);

}

FASTA::FASTA(FASTA&& other) noexcept
    : _name(other._name),
      _hts_fai(other._hts_fai),
      _size(other._size) {
    other._hts_fai = nullptr;
}

FASTA& FASTA::operator=(FASTA&& other) noexcept {

    if (this != &other) {
        _hts_fai = other._hts_fai;
        _size    = other._size;
        other._hts_fai = nullptr;
    }

    return *this;

}

FASTA::FASTA(const FASTA& other)
    : _name(other._name),
      _hts_fai(_open_fai(other._name)),
      _size(other._size) {}

FASTA::~FASTA() {

    _close_fai(_hts_fai);

}

faidx_t* FASTA::ptr() const {

    return _hts_fai;

}

std::string FASTA::name() const {

    return _name;

}

size_t FASTA::size() const {

    return _size;

}

hts_pos_t FASTA::length(hts_pos_t ix) const {

    _verify_index(*this, ix);
    return _get_length(*this, ix);

}

std::string FASTA::sequence(hts_pos_t ix) const {

    _verify_index(*this, ix);
    return _get_sequence(*this, ix);

}

seq_t FASTA::operator[](hts_pos_t ix) const {

    std::string _sequence = sequence(ix);
    return _as_int(_sequence);

}

void FASTA::as_int(hts_pos_t ix, view_t<base_t, DS_DIMS> buffer) const {

    std::string _sequence = sequence(ix);
    _as_int(_sequence, buffer);

}

hts_pos_t FASTA::max_length() const {

    hts_pos_t _max_length = 0;
    for (size_t ix = 0; ix < size(); ix++) {
        hts_pos_t len = length(ix);
        if (len > _max_length) {
            _max_length = len;
        }
    }

    return _max_length;

}

void FASTA::to_hdf5(HDF5::File& hdf5, const MPI::Manager& mpi) const {

    _fasta_to_hdf5<base_t>(*this, hdf5, mpi);

}



//
// CIGAR
//



static inline std::string _cigar_hts_format(CIGAR_t type) {

    switch (type) {
        case CIGAR_t::MATCH: {
            return "M";
        }
        case CIGAR_t::INS:   {
            return "I";
        }
        case CIGAR_t::DEL:   {
            return "D";
        }
        default:             {
            return "?";
        }
    }

}



CIGAR_op::CIGAR_op(CIGAR_t type) : _type(type), _length(1) {}

CIGAR_op::CIGAR_op(CIGAR_t type, hts_pos_t length) : _type(type), _length(length) {}

CIGAR_t CIGAR_op::type() const {

    return _type;

}

bool CIGAR_op::is_mod() const {

    return _is_mod(_type);

}

hts_pos_t CIGAR_op::length() const {

    return _length;

}

hts_pos_t CIGAR_op::rlength() const {

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

hts_pos_t CIGAR_op::qlength() const {

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

void CIGAR_op::extend(hts_pos_t val) {

    _length += val;

}

hts_pos_t CIGAR_op::pos() const {

    return _pos;

}

bool CIGAR_op::empty() const {

    return _length == 0;

}

bool CIGAR_op::advance(ssize_t offset) {

    if (_pos < _length + offset) {
        _pos++;
        return true;
    }
    return false;

}

std::pair<CIGAR_op, CIGAR_op> CIGAR_op::split(hts_pos_t n) const {

    hts_pos_t l1 = std::min(_length, n);
    hts_pos_t l2 = _length - l1;

    return std::pair<CIGAR_op, CIGAR_op>(
        CIGAR_op(_type, l1),
        CIGAR_op(_type, l2)
    );

}

CIGAR_op CIGAR_op::match() const {
    return CIGAR_op(CIGAR_t::MATCH, _length);
}

CIGAR::CIGAR(size_t size) {

    _str.reserve(size);

}

hts_pos_t CIGAR::rlength() const {

    hts_pos_t _rlength = 0;
    for (const auto& op : _str) {
        _rlength += op.rlength();
    }

    return _rlength;

}


std::string CIGAR::str() const {

    std::string _out_str;

    for (const auto& op: _str) {
        _out_str += op.str();
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

size_t CIGAR::size() const {

    return _str.size();

}

size_t CIGAR::bases() const {

    size_t bases = 0;
    for (const auto& op : _str) {
        bases += op.length();
    }

    return bases;

}

CIGAR_op CIGAR::operator[](size_t ix) const {

    return _str[ix];

}

auto CIGAR::begin() -> decltype(_str.begin()) {

    return _str.begin();

}
auto CIGAR::begin() const -> decltype(_str.begin()) {

    return _str.begin();

}
auto CIGAR::end() -> decltype(_str.end()) {

    return _str.end();

}
auto CIGAR::end() const -> decltype(_str.end()) {

    return _str.end();

}



//
// Alignment
//



Alignment::Alignment(
    htsFile* _hts_file,
    hts_itr_t* _hts_iter,
    const seq_t& _reference
) : _hts_file(_hts_file), _hts_iter(_hts_iter), _reference(_reference) {

    _hts_aln = _open_aln();

}

Alignment::Alignment(Alignment&& other) noexcept
    : _hts_file(other._hts_file),
      _hts_aln(other._hts_aln),
      _hts_iter(other._hts_iter),
      _reference(std::move(other._reference)) {

    other._hts_file = nullptr;
    other._hts_iter = nullptr;
    other._hts_aln  = nullptr;

}

Alignment& Alignment::operator=(Alignment&& other) noexcept {

    if (this != &other) {
        _hts_file = other._hts_file;
        _hts_iter = other._hts_iter;
        _hts_aln  = other._hts_aln;
        _reference = std::move(other._reference);
        other._hts_file = nullptr;
        other._hts_iter = nullptr;
        other._hts_aln  = nullptr;
    }

    return *this;

}

Alignment::~Alignment() {

    _close_aln(_hts_aln);
    _close_iter(_hts_iter);
    // _close_file(_hts_file);

}

bool Alignment::next() const {

    return sam_itr_next(
        _hts_file,
        _hts_iter,
        _hts_aln
    ) >= 0;

}

int64_t Alignment::reads() const {

    int64_t _reads = 0;
    while (next()) {
        if (aligned()) { _reads++; }
    }

    return _reads;

}

base_t Alignment::qbase(hts_pos_t ix) const {

    uint8_t* _query = bam_get_seq(_hts_aln);
    uint8_t _hts_base = bam_seqi(_query, ix);

    return _base_from_hts(_hts_base);

}

base_t Alignment::rbase(hts_pos_t ix) const {

    return _reference[ix];

}

hts_pos_t Alignment::qlength() const {

    return _hts_aln->core.l_qseq;

}

hts_pos_t Alignment::rlength() const {

    return _reference.size();

}

const std::vector<int8_t>& Alignment::reference() const {

    return _reference;

}

hts_pos_t Alignment::offset() const {

    return _hts_aln->core.pos;

}

uint8_t Alignment::mapq() const {

    return _hts_aln->core.qual;

}

std::span<const uint8_t> Alignment::phred() const {

    uint8_t* _quality = bam_get_qual(_hts_aln);
    return std::span<const uint8_t>(_quality, _quality + qlength());

}

template <typename dtype>
std::vector<dtype> Alignment::mask(uint8_t min, hts_pos_t window) const {

    return _get_mask<dtype>(phred(), min, window);

}

bool Alignment::aligned() const {

    return _hts_aln->core.tid >= 0;

}

CIGAR Alignment::cigar() const {

    std::string _md_tag = _get_md_tag(_hts_aln);
    std::span<const uint32_t> _cigar = _get_cigar_str(_hts_aln);

    return _get_cigar(_cigar, _md_tag);

}

bool Alignment::empty() const {

    return _hts_iter->finished;

}







//
// File
//



File::File(
    const std::string& filename,
    const FASTA& fasta,
    const MPI::Manager& mpi
) : _name(filename), _fasta(fasta), _mpi(mpi) {

    reset();

}

File::File(File&& other) noexcept
    : _hts_file(other._hts_file),
      _hts_index(other._hts_index),
      _hts_header(other._hts_header),
      _name(std::move(other._name)),
      _fasta(std::move(other._fasta)),
      _mpi(std::move(other._mpi))
{

    other._hts_file   = nullptr;
    other._hts_index  = nullptr;
    other._hts_header = nullptr;

}

File::~File() {

    _close_file(_hts_file);
    _close_index(_hts_index);
    _close_header(_hts_header);

}

const MPI::Manager& File::mpi() const {
    return _mpi;
}

htsFile* File::ptr() const {

    return _hts_file;

}

hts_idx_t* File::index() const {

    return _hts_index;

}

bam_hdr_t* File::header() const {

    return _hts_header;

}

htsFile* File::handle() const {

    _throw_if_not_exists(_name);
    return _open_file(_name);

}

hts_itr_t* File::iter(hts_pos_t ix) const {

    return _open_iter(*this, ix);

}

std::string File::name() const {

    return _name;

}

FileType File::type() const {

    return _filetype_from_hts(_hts_file->format.format);

}

const FASTA& File::fasta() const {

    return _fasta;

}

void File::reset() {

    _close_file(_hts_file);
    _close_header(_hts_header);
    _close_index(_hts_index);

    _throw_if_not_exists(_name);
    _sort(_name, _mpi);

    _hts_file   = _open_file(_name);
    _hts_header = _open_header(_hts_file, _name);
    _hts_index  = _open_index(*this);

    _throw_if_bad_lengths(*this);

    _set_reference(*this, _fasta.name());

}

std::string File::stem() const {

    return _stem(_name);

}

std::string File::path() const {

    return _path(_name);

}

seq_t File::reference(hts_pos_t ix) const {

    return _fasta[ix];

}

int64_t File::size() const {

    return hts_idx_nseq(_hts_index);

}

int64_t File::reads(hts_pos_t ix) const {

    if (type() == FileType::CRAM) {
        throw std::runtime_error("Read counts for individual references cannot be obtained for CRAM files.");
    } else {
        return _sam_bam_aligned(_hts_index, ix);
    }

}

int64_t File::reads() const {

    if (type() == FileType::CRAM) {
        return _cram_aligned(_hts_file, _hts_header);
    } else {
        return _sam_bam_aligned(_hts_index);
    }

}

int64_t File::unaligned_reads() const {

    return _unaligned(_hts_index);

}

hts_pos_t File::length(int64_t ix) const {

    return _hts_header->target_len[ix];

}

hts_pos_t File::max_length() const {

    hts_pos_t _max_length = 0;
    for (hts_pos_t ix = 0; ix < size(); ix++) {
        hts_pos_t len = length(ix);
        if (len > _max_length) {
            _max_length = len;
        }
    }

    return _max_length;

}

Alignment File::alignment(hts_pos_t ix) const {

    // htsFile* _aln_hts_file = _open_file(_name);
    hts_itr_t* _aln_hts_iter = iter(ix);
    seq_t _reference = _fasta[ix];

    return Alignment(_hts_file, _aln_hts_iter, _reference);

}



//
// FileGroup
//



FileGroup::FileGroup(
    const std::vector<std::string>& filenames,
    const HTS::FASTA& fasta,
    const MPI::Manager& mpi
) {

    _throw_if_has_duplicates(filenames);

    for (auto& name : filenames) {
        group.emplace_back(name, fasta, mpi);
    }

}

int64_t FileGroup::size() const {

    return group.size();

}

int64_t FileGroup::reads() {

    int64_t _reads = 0;
    for (auto& file : group) {
        _reads += file.reads();
    }

    return _reads;

}

int64_t FileGroup::unaligned_reads() const {

    int64_t _reads = 0;
    for (const auto& file : group) {
        _reads += file.unaligned_reads();
    }

    return _reads;

}

// All files are guaranteed to have the same number of references, and references of the same length, as compared to the FASTA.

int64_t FileGroup::references() const {

    if (size() > 0) {
        return group[0].size();
    } else {
        return 0;
    }

}

int64_t FileGroup::max_length() const {

    if (size() > 0) {
        return group[0].max_length();
    } else {
        return 0;
    }

}

auto FileGroup::begin() -> decltype(group.begin()) {

    return group.begin();

}

auto FileGroup::begin() const -> decltype(group.begin()) {

    return group.begin();

}

auto FileGroup::end() -> decltype(group.end()) {

    return group.end();

}

auto FileGroup::end() const -> decltype(group.end()) {

    return group.end();

}



//
// Template specialisation
//



template std::vector<float> Alignment::mask<float>(uint8_t, hts_pos_t) const;



} // namespace HTS
