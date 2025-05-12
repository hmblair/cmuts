#include "hts.hpp"

namespace HTS {

void __disable_logging() {

    hts_set_log_level(htsLogLevel::HTS_LOG_OFF);

}

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

static inline seq_t _to_int(
    const std::string& sequence
) {
    int64_t length = sequence.length();
    seq_t result(length);
    for (int64_t i = 0; i < length; i++) {
        result[i] = _to_int(sequence[i]);
    }
    return result;
}

static inline void _to_int(
    const std::string& sequence,
    view_t<base_t, DS_DIMS> buffer
) {
    for (int64_t i = 0; i < sequence.length(); i++) {
        buffer(i) = _to_int(sequence[i]);
    }
}

static inline size_t _bytes_from_seq(size_t n) {

    return (n + 3) / 4;

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



//
// Simple FASTA tools
//




static inline void _write_fasta_record(const std::string& name, const std::string& sequence, std::ofstream& file) {

    file << ">" << name << "\n";
    file << sequence << "\n";

}

FASTA::FASTA(const std::string& filename) : _name(filename), _out_file(filename) {}

void FASTA::write(const std::string& name, const std::string& sequence) {

    _write_fasta_record(name, sequence, _out_file);

}



//
// Binary FASTA tools
//



static inline std::vector<uint8_t> _to_binary(const std::string& sequence) {

    size_t len = sequence.length();
    size_t bytes = _bytes_from_seq(len);
    std::vector<uint8_t> binary(bytes);

    for (size_t i = 0; i < len; i++) {

        base_t nuc = _to_int(sequence[i]);
        if (nuc == IX_UNK) {
            throw std::runtime_error("Invalid nucleotide \"" + std::to_string(sequence[i]) + "\" encountered.");
        }

        binary[i / 4] |= (nuc << ((3 - (i % 4)) * 2));

    }

    return binary;

}

static inline uint8_t _from_binary(uint8_t binary, hts_pos_t j) {

    return (binary >> ((3 - j) * 2)) & 0x3;

}

hts_pos_t _mod_round_up(hts_pos_t val, hts_pos_t div) {

    hts_pos_t out = val % div;
    if (out == 0) { out = div; }
    return out;

}

static inline seq_t _from_binary(const std::vector<uint8_t>& binary, size_t length) {

    seq_t sequence(length);

    hts_pos_t ix, jx;
    for (ix = 0; ix < binary.size() - 1; ix++) {
        for (jx = 0; jx < 4; jx++) {
            sequence[4 * ix + jx] = _from_binary(binary[ix], jx);
        }
    }

    hts_pos_t end = _mod_round_up(length, 4);
    for (jx = 0; jx < end; jx++) {
        sequence[4 * ix + jx] = _from_binary(binary[ix], jx);
    }

    return sequence;

}

static inline bool _is_fasta_header(const std::string& line) {

    return !line.empty() && line[0] == '>';

}

static inline std::string _get_fasta_sequence(std::ifstream& file) {

    std::string line, sequence;
    while (std::getline(file, line) && !_is_fasta_header(line)) {
        sequence += line;
    }

    return sequence;

}




static inline HeaderUnit _read_header_unit(std::ifstream& file) {

    HeaderUnit unit;
    file.read(reinterpret_cast<char*>(&unit.length), sizeof(hts_pos_t));
    file.read(reinterpret_cast<char*>(&unit.sequences), sizeof(hts_pos_t));
    return unit;

}

static inline std::vector<HeaderUnit> _read_header(const std::string& name) {

    std::ifstream file(name);
    std::vector<HeaderUnit> units;
    HeaderUnit unit;

    while (unit.length != EOH) {
        unit = _read_header_unit(file);
        if (unit.length > 0) {
            units.push_back(unit);
        }
    }

    file.close();
    return units;

}


static inline void _write_header_unit(std::ofstream& file, const HeaderUnit& unit) {

    if (unit.sequences > 0) {
        file.write(reinterpret_cast<const char*>(&unit.length), sizeof(hts_pos_t));
        file.write(reinterpret_cast<const char*>(&unit.sequences), sizeof(hts_pos_t));
    }

}

static inline void _write_eoh(std::ofstream& file) {

    file.write(reinterpret_cast<const char*>(&EOH), sizeof(hts_pos_t));

}

static inline void _write_header(const std::string& name, const std::vector<HeaderUnit>& units) {

    // _throw_if_exists(name);
    std::ofstream file(name);

    for (const auto& unit : units) {
        _write_header_unit(file, unit);
    }
    _write_eoh(file);

    file.close();

}

static inline void _write_sequence(std::ofstream& file, const std::string& sequence) {

    std::vector<uint8_t> binary = _to_binary(sequence);
    file.write(reinterpret_cast<const char*>(binary.data()), binary.size() * sizeof(uint8_t));

}

static inline void _write_sequences(const std::string& fasta, const std::string& out) {

    if (fasta == out) {
        throw std::runtime_error("The input and output files must be different.");
    }
    _throw_if_not_exists(fasta);
    // _throw_if_exists(out);

    std::ifstream ifile(fasta);
    std::ofstream ofile(out, std::ios::binary | std::ios::app);

    std::string line;
    std::getline(ifile, line);
    std::string sequence = _get_fasta_sequence(ifile);

    if (sequence.empty()) { return; }
    _write_sequence(ofile, sequence);

    while (!sequence.empty()) {
        sequence = _get_fasta_sequence(ifile);
        _write_sequence(ofile, sequence);
    }

    ifile.close();
    ofile.close();

}

static inline seq_t _read_sequence(std::ifstream& file, size_t length, size_t offset) {

    size_t bytes = _bytes_from_seq(length);
    std::vector<uint8_t> binary(bytes);

    file.seekg(offset);
    file.read(reinterpret_cast<char*>(binary.data()), bytes * sizeof(uint8_t));

    return _from_binary(binary, length);

}

static inline std::vector<HeaderUnit> _get_header(const std::string& name) {

    _throw_if_not_exists(name);
    std::ifstream file(name, std::ios::in);
    std::vector<HeaderUnit> units;
    HeaderUnit unit;

    std::string line;
    std::getline(file, line);
    std::string sequence = _get_fasta_sequence(file);
    if (sequence.empty()) { return units; }

    size_t length = sequence.length();
    unit.length = length;
    unit.sequences = 1;

    while (!sequence.empty()) {

        sequence = _get_fasta_sequence(file);
        size_t length = sequence.length();
        if (length == unit.length) {
            unit.sequences++;
        } else {
            units.push_back(unit);
            unit.length = length;
            unit.sequences = 1;
        }

    }

    if (unit.sequences > 0 && unit.length > 0) {
        units.push_back(unit);
    }

    return units;

}

template <typename dtype>
static inline HDF5::Memspace<dtype, DS_DIMS> __memspace(
    const HTS::BinaryFASTA& fasta,
    HDF5::File& hdf5
) {
    std::vector<size_t> dims = {fasta.size(), static_cast<size_t>(fasta.max_length())};
    return hdf5.memspace<dtype, DS_DIMS>(dims, SEQUENCE_DS);
}

template <typename dtype>
static inline void _fasta_to_hdf5(
    HTS::BinaryFASTA& fasta,
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
                fasta.put(ix + jx, arr);
            }
        }
        memspace.safe_write(ix);
        memspace.clear();
    }

}



static inline void _generate_binary_index(const std::string& name) {

    std::string _name = _path(name) + CMUTS_FASTA;
    std::ofstream file(_name);
    std::vector<HeaderUnit> units = _get_header(name);
    _write_header(_name, units);
    _write_sequences(name, _name);
    file.close();

}




BinaryFASTA::BinaryFASTA(const std::string& name, const MPI::Manager& mpi) {

    if (mpi.root()) { _generate_binary_index(name); }
    mpi.barrier();
    _name = _path(name) + CMUTS_FASTA;

    units   = _read_header(_name);
    _offset = (units.size() * 2 + 1) * sizeof(hts_pos_t);
    file = std::ifstream(_name);

}

BinaryFASTA::BinaryFASTA(BinaryFASTA&& other) noexcept
    : _name(std::move(other._name)),
      units(std::move(other.units)),
      file(std::move(other.file)),
      _offset(other._offset) {}

BinaryFASTA& BinaryFASTA::operator=(BinaryFASTA&& other) noexcept {

    if (this != &other) {
        _name = other._name;
        units = other.units;
        file = std::move(other.file);
        _offset = other._offset;
    }

    return *this;

}

BinaryFASTA::BinaryFASTA(const BinaryFASTA& other)
    : _name(other._name),
      units(other.units),
      file(other._name),
      _offset(other._offset) {}

std::string BinaryFASTA::name() const {

    return _name;

}

size_t BinaryFASTA::length(hts_pos_t ix) const {

    hts_pos_t count = 0;
    for (const auto& unit : units) {
        if (count + unit.sequences > ix) {
            return unit.length;
        }
        count += unit.sequences;
    }

    return -1;

}

size_t BinaryFASTA::size() const {

    size_t _sequences = 0;
    for (const auto& unit : units) {
        _sequences += unit.sequences;
    }

    return _sequences;

}

std::pair<size_t, size_t> BinaryFASTA::offset(hts_pos_t ix) const {

    size_t val = 0;
    size_t length;
    hts_pos_t count = 0;
    for (const auto& unit : units) {

        size_t bytes = _bytes_from_seq(unit.length);

        if (count + unit.sequences > ix) {
            val += static_cast<size_t>(ix - count) * bytes;
            length = unit.length;
            break;
        } else {
            val += unit.sequences * bytes;
            count += unit.sequences;
        }

    }

    // Account for the header
    size_t __offset = val * sizeof(uint8_t) + _offset;
    return std::pair<size_t, size_t>(__offset, length);

}

seq_t BinaryFASTA::operator[](hts_pos_t ix) {

    auto [val, length] = offset(ix);
    return _read_sequence(file, length, val);

}

void BinaryFASTA::put(hts_pos_t ix, view_t<base_t, DS_DIMS> buffer) {

    seq_t _seq = operator[](ix);
    for (hts_pos_t jx = 0; jx < _seq.size(); jx++) {
        buffer(jx) = _seq[jx];
    }

}

hts_pos_t BinaryFASTA::max_length() const {

    hts_pos_t _max_length = 0;
    for (size_t ix = 0; ix < size(); ix++) {
        hts_pos_t len = length(ix);
        if (len > _max_length) {
            _max_length = len;
        }
    }

    return _max_length;

}

void BinaryFASTA::to_hdf5(HDF5::File& hdf5, const MPI::Manager& mpi) {

    _fasta_to_hdf5<base_t>(*this, hdf5, mpi);

}







// Other







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


static inline void _close_bgzf(BGZF*& _bgzf_file) {

    if (_bgzf_file != nullptr) {
        bgzf_close(_bgzf_file);
        _bgzf_file = nullptr;
    }

}

static inline void _validate_bgzf(BGZF* _bgzf_file) {

    char magic[4];
    if (bgzf_read(_bgzf_file, magic, 4) != 4 || strncmp(magic, "BAM\1", 4) != 0) {
        _close_bgzf(_bgzf_file);
        throw std::runtime_error("Magic check failed -- invalid BAM file.");
    }

}

static inline BGZF* _open_bgzf(const std::string& name) {

    BGZF* _bgzf_file = bgzf_open(name.c_str(), "r");
    if (_bgzf_file == nullptr) {
        throw std::runtime_error("Failed to open the file \"" + name + "\" using BGZF.");
    }

    _validate_bgzf(_bgzf_file);
    return _bgzf_file;

}

static inline hts_pos_t _get_header_length(BGZF* _bgzf_file) {

    int32_t length;
    if (bgzf_read(_bgzf_file, &length, sizeof(length)) != sizeof(length)) {
        _close_bgzf(_bgzf_file);
        throw std::runtime_error("Failed to get the length of the header.");
    }

    return static_cast<hts_pos_t>(length);

}

static inline std::string _read_bgzf_line(
    BGZF* _bgzf_file,
    char* buffer,
    hts_pos_t& ix,
    hts_pos_t length
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

static inline bool _is_sorted_rec(const std::string& str) {

    size_t pos = str.find(HEADER_SORT_KEY) + HEADER_SORT_KEY.length();
    std::string value = str.substr(pos);

    if (value == BGZF_SORTED) {
        return true;
    }
    return false;

}

static inline bool _is_sorted(const std::string& name) {

    BGZF* _bgzf_file = _open_bgzf(name);

    char buffer[BGZF_BUFFER];
    hts_pos_t ix = 0;
    hts_pos_t length = _get_header_length(_bgzf_file);
    std::string line = _read_bgzf_line(_bgzf_file, buffer, ix, length);

    _close_bgzf(_bgzf_file);
    return _is_sorted_rec(line);

}

static inline hts_pos_t _length_from_rec(const std::string& str) {

    size_t pos = str.find(HEADER_LEN_KEY) + HEADER_LEN_KEY.length();
    std::string value = str.substr(pos);
    return static_cast<hts_pos_t>(std::stoi(value));

}

static inline bool _is_seq_rec(const std::string& line) {

    return line.starts_with("@SQ");

}

static inline std::vector<hts_pos_t> _get_reference_lengths(const std::string& name) {

    BGZF* _bgzf_file = _open_bgzf(name);
    hts_pos_t length = _get_header_length(_bgzf_file);

    std::vector<hts_pos_t> sizes;

    // Read header line by line
    char buffer[BGZF_BUFFER];
    hts_pos_t ix = 0;
    while (ix < length) {
        std::string line = _read_bgzf_line(_bgzf_file, buffer, ix, length);
        if (_is_seq_rec(line)) {
            sizes.push_back(_length_from_rec(line));
        }
    }

    _close_bgzf(_bgzf_file);
    return sizes;

}

static inline void _throw_if_not_sorted(const std::string& name) {

    if (!_is_sorted(name)) {
        throw std::runtime_error("The input file \"" + name + "\" is not sorted.");
    }

}

static inline void _sort(const std::string& name, const MPI::Manager& mpi) {

    mpi.barrier();

    if (mpi.root() && !_is_sorted(name)) {

        mpi.out() << "        Sorting \"" << name << "\".\n";

        std::string unsorted_name = _path(name) + "_UNSORTED.bam";
        _safe_move(name, unsorted_name);

        std::string command = "samtools sort -o " + name + " -@ " + std::to_string(mpi.size()) + " " + unsorted_name + " > cmuts.log 2>&1";
        if (system(command.c_str()) != 0) {
            _safe_move(unsorted_name, name);
            throw std::runtime_error("There was an error sorting the file \"" + name + "\".");
        }
        mpi.divide();

    }

    mpi.barrier();

}

static inline hts_idx_t* _open_index(const File& file) {

    _throw_if_not_sorted(file.name());

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

void _throw_if_bad_lengths(File& file) {

    BinaryFASTA& fasta = file.fasta();
    if (file.size() != fasta.size()) {
        throw std::runtime_error("The number of references in the header for the file \"" + file.name() + "\" (" + std::to_string(file.size()) + ") does not match the number of sequences in the FASTA (" + std::to_string(fasta.size()) + ").");
    }

    std::vector<hts_pos_t> _lengths = file.lengths();
    for (hts_pos_t ix = 0; ix < file.size(); ix++) {
        if (_lengths[ix] != fasta.length(ix)) {
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



//
// CIGAR
//



static inline std::string _cigar_hts_format(CIGAR_t type) {

    switch (type) {
        case CIGAR_t::MATCH: { return "M"; }
        case CIGAR_t::INS:   { return "I"; }
        case CIGAR_t::DEL:   { return "D"; }
        default:             { return "?"; }
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

CIGAR_op CIGAR::back() const {

    if (!_str.empty()) { return _str.back(); }
    return CIGAR_op();

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
    const BinaryFASTA& fasta,
    const MPI::Manager& mpi
) : _name(filename), _fasta(fasta), _mpi(mpi) {

    reset();

}

File::File(File&& other) noexcept
    : _hts_file(other._hts_file),
      _hts_index(other._hts_index),
      _name(std::move(other._name)),
      _fasta(std::move(other._fasta)),
      _mpi(std::move(other._mpi))
{

    other._hts_file   = nullptr;
    other._hts_index  = nullptr;

}

File::~File() {

    _close_file(_hts_file);
    _close_index(_hts_index);

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

BinaryFASTA& File::fasta() {

    return _fasta;

}

void File::reset() {

    _close_file(_hts_file);
    _close_index(_hts_index);

    _throw_if_not_exists(_name);
    _sort(_name, _mpi);

    _hts_file   = _open_file(_name);
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

seq_t File::reference(hts_pos_t ix) {

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
        throw std::runtime_error("Not implemented yet.");
        // return _cram_aligned(_hts_file, _hts_header);
    } else {
        return _sam_bam_aligned(_hts_index);
    }

}

int64_t File::unaligned_reads() const {

    return _unaligned(_hts_index);

}

std::vector<hts_pos_t> File::lengths() const {

    return _get_reference_lengths(_name);

}

hts_pos_t File::max_length() const {

    std::vector<hts_pos_t> _lengths = lengths();
    return *std::max_element(_lengths.begin(), _lengths.end());

}

Alignment File::alignment(hts_pos_t ix) {

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
    const HTS::BinaryFASTA& fasta,
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
