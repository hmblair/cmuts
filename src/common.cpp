#include "common.hpp"
#include <stdexcept>
#include <string>





//
// General file utilities
//





bool _exists(const std::string& filename) {

    return std::filesystem::exists(filename);

}



void _throw_if_exists(const std::string& filename) {

    if (_exists(filename)) {
        throw std::runtime_error("The file \"" + filename + "\" already exists.");
    }

}


void _throw_if_not_exists(const std::string& filename) {

    if (!_exists(filename)) {
        throw std::runtime_error("The file \"" + filename + "\" does not exist.");
    }

}


void _safe_move(const std::string& src, const std::string& dst) {

    _throw_if_not_exists(src);
    _throw_if_exists(dst);
    std::filesystem::rename(src, dst);

}


std::string _stem(const std::string& name) {

    std::filesystem::path path(name);
    return path.stem().string();

}


std::string _path(const std::string& name) {

    std::filesystem::path path(name);
    std::string stem = path.stem().string();
    std::string parent = path.parent_path().string();
    if (parent.empty()) { return stem; }
    return parent + "/" + stem;

}


bool _has_duplicate_paths(const std::vector<std::string>& paths) {

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


void _throw_if_has_duplicate_paths(const std::vector<std::string>& paths) {

    if (_has_duplicate_paths(paths)) {
        throw std::runtime_error("Duplicate input paths (modulo extensions) are not allowed.");
    }

}





//
// BGZF file utilities
//





BGZF* _open_bgzf(const std::string& name) {

    _throw_if_not_exists(name);

    BGZF* _bgzf_file = bgzf_open(name.c_str(), "r");
    if (_bgzf_file == nullptr) {
        throw std::runtime_error("Failed to open the file \"" + name + "\" using BGZF.");
    }

    return _bgzf_file;

}


void _close_bgzf(BGZF*& _bgzf_file) {

    if (_bgzf_file != nullptr) {
        bgzf_close(_bgzf_file);
        _bgzf_file = nullptr;
    }

}


void _read_bgzf(BGZF* _bgzf_file, void* buffer, int32_t size) {

    if (_bgzf_file == nullptr) {
        throw std::runtime_error("Null BGZF pointer.");
    }

    auto bytes = static_cast<int32_t>(
        bgzf_read(_bgzf_file, buffer, size)
    );

    if (bytes != size) {

        _close_bgzf(_bgzf_file);
        throw std::runtime_error("BGZF failed to read " + std::to_string(size - bytes) + " of " + std::to_string(size) + " bytes.");

    }

}


std::vector<uint8_t> _read_bgzf(BGZF* _bgzf_file, int32_t size) {

    std::vector<uint8_t> buffer(size);
    _read_bgzf(_bgzf_file, buffer.data(), size);
    return buffer;

}


void _seek_bgzf(BGZF* _bgzf_file, int64_t ptr) {

    if (bgzf_seek(_bgzf_file, ptr, SEEK_SET) < 0) {

        _close_bgzf(_bgzf_file);
        throw std::runtime_error("Failed to seek to position " + std::to_string(ptr) + ".");

    }

}


template <typename dtype>
dtype _read_bgzf_line(BGZF* _bgzf_file, uint8_t delimiter) {

    kstring_t tmp = KS_INITIALIZE;

    int length = bgzf_getline(_bgzf_file, static_cast<int>(delimiter), &tmp);
    if (length < -1) {
        throw std::runtime_error("Failure in _read_bgzf_line.");
    }

    dtype out(tmp.s, tmp.s + length);
    ks_free(&tmp);
    return out;

}


template <typename dtype>
dtype _read_bgzf_single(BGZF* _bgzf_file) {

    dtype value;
    _read_bgzf(_bgzf_file, &value, sizeof(dtype));
    return value;

}





//
// zlib file utilities
//





static inline std::string _get_stream_msg(z_stream* stream) {

    if (stream->msg == nullptr) {
        return "(no information given)";
    } else {
        return {stream->msg, strlen(stream->msg)};
    }

}


static inline z_stream _open_z_stream(std::span<const uint8_t> data) {

    z_stream stream;

    stream.zalloc    = Z_NULL;
    stream.zfree     = Z_NULL;
    stream.opaque    = Z_NULL;
    stream.next_in   = const_cast<uint8_t*>(data.data());
    stream.avail_in  = data.size();
    stream.total_in  = 0;
    stream.next_out  = nullptr;
    stream.avail_out = 0;
    stream.total_out = 0;

    int err = inflateInit2(&stream, ZLIB_WINDOW);
    if (err != Z_OK) {
        throw std::runtime_error("Failed to initialize z_stream: " + _get_stream_msg(&stream));
    }

    return stream;

}


static inline void _close_z_stream(z_stream& stream) {

    (void)inflateEnd(&stream);

}


static inline void _read_z_stream(z_stream* stream, std::vector<uint8_t>& buffer) {

    if (stream == nullptr) {
        throw std::runtime_error("Null pointer in _read_z_stream.");
    }

    stream->total_out = 0;
    stream->total_in  = 0;

    stream->next_out  = buffer.data();
    stream->avail_out = buffer.size();

    int err = inflate(stream, Z_NO_FLUSH);
    if (err != Z_OK && err != Z_STREAM_END) {

        throw std::runtime_error(
            "Failed to read from the zlib stream: " + _get_stream_msg(stream)
        );

    }

}





//
// FASTA
//





static inline bool _is_fasta_header(const std::string& line) {

    return !line.empty() && line[0] == '>';

}


static inline std::string _get_fasta_sequence(std::fstream& file) {

    std::string line;
    std::string sequence;

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





//
// ByteStream
//





template <typename dtype>
static inline int32_t _find_first_occurrence(
    const dtype* begin,
    const dtype* end,
    dtype value
) {

    auto it = std::find(begin, end, value);
    if (it == end) {
        return -1;
    } else {
        return std::distance(begin, it);
    }

}





ByteStream::ByteStream(int32_t size) : _size(size), _remaining(size) {}


std::string ByteStream::str(uint8_t delimiter) {

    std::vector<uint8_t> values = line(delimiter);
    return {values.begin(), values.end()};

}


int32_t ByteStream::size() const {

    return _size;

}


int32_t ByteStream::remaining() const {

    return _remaining;

}


bool ByteStream::end() const noexcept {

    return (_remaining <= 0);

}


int32_t ByteStream::bits(int32_t length) {

    int32_t out            = 0;
    int32_t bits_remaining = length;

    while (bits_remaining > 0) {

        if (_n_bits == 0) { _byte = byte(); _n_bits = BYTE; }

        int32_t _to_read = std::min(bits_remaining, static_cast<int32_t>(_n_bits));

        // Make room for the bits we are going to read

        out <<= _to_read;

        // Read the desired number of bits into out

        uint8_t tmp = (_byte & BIT_MASK[_to_read]);
        out |= (tmp >> (BYTE - _to_read));

        // Update accordingly

        _byte   <<= _to_read;
        _n_bits -= _to_read;
        bits_remaining -= _to_read;

    }

    return out;

}


void ByteStream::memcpy(uint8_t* dest, int32_t n) {

    std::vector<uint8_t> values = bytes(n);
    std::memcpy(dest, values.data(), n);

}


int32_t ByteStream::itf8() {

    // The first byte tells us the number of remaining bytes

    uint8_t first = byte();
    int32_t size = std::countl_one(static_cast<uint8_t>(first & 0xF0));

    // Initialize the value with the least 7 - size significant bits

    auto value = static_cast<int32_t>(first & ITF8_MASK[size]);
    if (size == 0) { return value; }

    // Read the appropriate number of bytes

    std::vector<uint8_t> remaining = bytes(size);

    // Compute the value by appending the data from the read bytes

    for (int32_t ix = 0; ix < size - 1; ix++) {
        value = (value << 8) | remaining[ix];
    }

    if (size == 4) {
        value = (value << 4) | (remaining[3] & 0x0F);
    } else {
        value = (value << 8) | remaining[size - 1];
    }

    return value;

}


int64_t ByteStream::ltf8() {

    // The first byte tells us the number of remaining bytes

    uint8_t first = byte();
    int32_t size = std::countl_one(first);

    // Initialize the value with all bits

    auto value = static_cast<int64_t>(first);
    if (size == 0) { return value; }

    // Read the appropriate number of bytes

    std::vector<uint8_t> remaining = bytes(size);

    // Compute the value by appending the data from the read bytes

    for (int32_t ix = 0; ix < size; ix++) {
        value = (value << 8) | remaining[ix];
    }

    return value;

}





//
// dataStream
//





dataStream::dataStream(std::span<const uint8_t> data)
    : ByteStream(static_cast<int32_t>(data.size())), _data(data) {}


uint8_t dataStream::byte() {

    uint8_t value = _data[_pos];
    _pos++;
    _remaining--;
    return value;

}


std::vector<uint8_t> dataStream::bytes(int32_t length) {

    if (length <= 0 || length + _pos > _data.size()) {
        throw std::runtime_error("Failed to get " + std::to_string(length) + " bytes from the dataStream.");
    }

    std::vector<uint8_t> values(length);
    std::memcpy(values.data(), _data.data() + _pos, length);

    _pos += length;
    _remaining -= length;

    return values;

}


std::vector<uint8_t> dataStream::line(uint8_t delimiter) {

    std::vector<uint8_t> line;
    line.reserve(LINE_RESERVE);

    int32_t pos = _find_first_occurrence<uint8_t>(
        _data.data() + _pos, _data.data() + _data.size(), delimiter
    );

    if (pos == -1) {
        line.insert(line.end(), _data.begin() + _pos, _data.end());
    } else if (pos > 0) {
        line.insert(line.end(), _data.begin() + _pos, _data.begin() + _pos + pos);
    }

    _pos += (pos + 1);
    _remaining  -= (pos + 1);

    return line;

}


std::string dataStream::str(uint8_t delimiter) {

    std::string line;
    line.reserve(LINE_RESERVE);

    int32_t pos = _find_first_occurrence<uint8_t>(
        _data.data() + _pos, _data.data() + _data.size(), delimiter
    );

    if (pos == -1) {
        line.insert(line.end(), _data.begin() + _pos, _data.end());
    } else if (pos > 0) {
        line.insert(line.end(), _data.begin() + _pos, _data.begin() + _pos + pos);
    }

    _pos += (pos + 1);
    _remaining  -= (pos + 1);

    return line;

}





//
// ransStream
//





static inline std::vector<uint8_t> _decompress_rans4x8(std::span<const uint8_t> data) {

    uint32_t _out_size = 0;
    uint8_t* _out = rans_uncompress(const_cast<uint8_t*>(data.data()), data.size(), &_out_size);
    if (_out == nullptr) {
        throw std::runtime_error("Failed to decompress the RANS4x8 stream.");
    }

    std::vector<uint8_t> buffer(_out, _out + _out_size);
    free(_out);

    return buffer;

}


ransStream::ransStream(std::vector<uint8_t> data)
    : dataStream(data), _rans_data(std::move(data)) {}


ransStream::ransStream(std::span<const uint8_t> data, int32_t raw)
    : ransStream(_decompress_rans4x8(data)) {

    if (_rans_data.size() != raw) {
        throw std::runtime_error("RANS4x8 decompresssion failed.");
    }

}





//
// zlibStream
//





zlibStream::zlibStream(std::span<const uint8_t> data, int32_t raw, int32_t buffer)
    : ByteStream(raw),
    _stream(_open_z_stream(data)),
    _buffer(buffer),
    _buffer_pos(_buffer.size()),
    _buffer_end(_buffer.size()),
    _open(true) {}


zlibStream::~zlibStream() {

    if (_open) {
        _close_z_stream(_stream);
        _open = false;
    }

}


void zlibStream::fill() {

    if (!_open) {
        throw std::runtime_error("The zlib stream is not open.");
    }

    _read_z_stream(&_stream, _buffer);
    _buffer_end = static_cast<int32_t>(_stream.total_out);
    _buffer_pos = 0;

}


uint8_t zlibStream::byte() {

    if (_buffer_pos >= _buffer_end) { fill(); }

    uint8_t value = _buffer[_buffer_pos];

    _remaining--;
    _buffer_pos++;

    return value;

}


std::vector<uint8_t> zlibStream::bytes(int32_t length) {

    std::vector<uint8_t> values(length);
    int32_t read = 0;

    while (read < length) {

        if (_buffer_pos >= _buffer_end) { fill(); }

        int32_t _to_read = std::min(length - read, static_cast<int32_t>(_buffer_end - _buffer_pos));
        if (_to_read <= 0) {
            throw std::runtime_error("Invalid number of bytes (" + std::to_string(_to_read) + ") in zlibStream.");
        }

        std::memcpy(values.data() + read, _buffer.data() + _buffer_pos, _to_read);

        _remaining  -= _to_read;
        _buffer_pos += _to_read;
        read        += _to_read;

    }

    return values;


}


std::vector<uint8_t> zlibStream::line(uint8_t delimiter) {

    if (_buffer_pos >= _buffer_end) { fill(); }

    std::vector<uint8_t> line;
    line.reserve(LINE_RESERVE);
    int32_t pos = _find_first_occurrence<uint8_t>(
        _buffer.data() + _buffer_pos, _buffer.data() + _buffer_end, delimiter
    );

    while (pos == -1 && !end()) {

        line.insert(
            line.end(),
            _buffer.begin() + _buffer_pos,
            _buffer.begin() + _buffer_end
        );

        _remaining -= (_buffer_end - _buffer_pos);

        fill();
        pos = _find_first_occurrence<uint8_t>(
            _buffer.data() + _buffer_pos, _buffer.data() + _buffer_end, delimiter
        );

    }

    if (pos == -1) {
        line.insert(
            line.end(),
            _buffer.begin() + _buffer_pos,
            _buffer.begin() + _buffer_end
        );
    } else if (pos > 0) {
        line.insert(
            line.end(),
            _buffer.begin() + _buffer_pos,
            _buffer.begin() + _buffer_pos + pos
        );
    }

    _buffer_pos += (pos + 1);
    _remaining  -= (pos + 1);

    return line;

}


std::string zlibStream::str(uint8_t delimiter) {

    if (_buffer_pos >= _buffer_end) { fill(); }

    std::string line;
    line.reserve(LINE_RESERVE);
    int32_t pos = _find_first_occurrence<uint8_t>(
        _buffer.data() + _buffer_pos, _buffer.data() + _buffer_end, delimiter
    );

    while (pos == -1 && !end()) {

        line.insert(
            line.end(),
            _buffer.begin() + _buffer_pos,
            _buffer.begin() + _buffer_end
        );

        _remaining -= (_buffer_end - _buffer_pos);

        fill();
        pos = _find_first_occurrence<uint8_t>(
            _buffer.data() + _buffer_pos, _buffer.data() + _buffer_end, delimiter
        );

    }

    if (pos == -1) {
        line.insert(
            line.end(),
            _buffer.begin() + _buffer_pos,
            _buffer.begin() + _buffer_end
        );
    } else if (pos > 0) {
        line.insert(
            line.end(),
            _buffer.begin() + _buffer_pos,
            _buffer.begin() + _buffer_pos + pos
        );
    }

    _buffer_pos += (pos + 1);
    _remaining  -= (pos + 1);

    return line;

}


void zlibStream::update(std::span<const uint8_t> data) {

    _stream.next_in  = const_cast<uint8_t*>(data.data());
    _stream.avail_in = data.size();

}





//
// bgzfFileStream
//





bgzfFileStream::bgzfFileStream(BGZF* file, int32_t size)
    : ByteStream(size), _file(file) {}


uint8_t bgzfFileStream::byte() {

    _remaining--;
    return _read_bgzf_single<uint8_t>(_file);

}


std::vector<uint8_t> bgzfFileStream::bytes(int32_t length) {

    _remaining -= length;
    return _read_bgzf(_file, length);

}


std::vector<uint8_t> bgzfFileStream::line(uint8_t delimiter) {

    auto line = _read_bgzf_line<std::vector<uint8_t>>(_file, delimiter);
    _remaining -= static_cast<int32_t>(line.size());
    return line;

}


std::string bgzfFileStream::str(uint8_t delimiter) {

    auto line = _read_bgzf_line<std::string>(_file, delimiter);
    _remaining -= static_cast<int32_t>(line.size());
    return line;

}




//
// zlibFileStream
//





zlibFileStream::zlibFileStream(BGZF* file, int32_t size, int32_t raw, int32_t buffer)
    : zlibStream(_buffer_span, raw, buffer),
      _bgzf(file, size),
      _bgzf_buffer(buffer),
      _bgzf_buffer_pos(buffer),
      _in_remaining(size) {}


void zlibFileStream::fill() {

    if (_bgzf_buffer_pos >= _bgzf_buffer.size()) {

        int32_t _to_read = std::min(_in_remaining, static_cast<int32_t>(_bgzf_buffer.size()));
        _bgzf_buffer = _bgzf.bytes(_to_read);
        _buffer_span = std::span<uint8_t>(_bgzf_buffer.begin(), _bgzf_buffer.end());

        update(_buffer_span);
        _in_remaining -= _to_read;
        _bgzf_buffer_pos = 0;

    }

    zlibStream::fill();
    _bgzf_buffer_pos += static_cast<int32_t>(_stream.total_in);

}





//
// Template specializations
//


template uint8_t _read_bgzf_single<uint8_t>(BGZF* file);
template int32_t _read_bgzf_single<int32_t>(BGZF* file);






namespace HTS {





void _disable_logging() {

    hts_set_log_level(htsLogLevel::HTS_LOG_OFF);

}


base_t _from_char(char base) {

    switch (base) {

        case A:  { return IX_A; }
        case C:  { return IX_C; }
        case G:  { return IX_G; }
        case T:  { return IX_T; }
        case U:  { return IX_U; }
        default: { return IX_UNK; }

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


std::string str(const seq_t& sequence) {

    std::string str(sequence.size(), '\0');
    for (hts_pos_t ix = 0; ix < sequence.size(); ix++) {
        str[ix] = _to_char(sequence[ix]);
    }

    return str;

}






//
// Header-related utilities
//





FileType _get_filetype(BGZF* _bgzf_file) {

    std::array<char, MAGIC_SIZE> magic = {0};
    _read_bgzf(_bgzf_file, magic.data(), magic.size());

    if (strncmp(magic.data(), SAM_MAGIC.c_str(), magic.size()) == 0) {

        return FileType::SAM;

    } else if (strncmp(magic.data(), BAM_MAGIC.c_str(), magic.size()) == 0) {

        return FileType::BAM;

    } else if (strncmp(magic.data(), CRAM_MAGIC.c_str(), magic.size()) == 0) {

        return FileType::CRAM;

    } else {

        _close_bgzf(_bgzf_file);
        throw std::runtime_error("Magic check failed -- invalid file.");

    }

}


static inline bool _is_sorted(const std::string& str) {

    int32_t pos = str.find(HEADER_SORT_KEY) + HEADER_SORT_KEY.length();
    std::string value = str.substr(pos);

    return (value == BGZF_SORTED);

}


static inline bool _is_seq(const std::string& line) {

    return line.starts_with("@SQ");

}

static inline bool _is_aux(const std::string& line) {

    return line.starts_with("@PG");

}


Header _read_sam_header(std::unique_ptr<ByteStream>& block, int32_t length) {

    Header header;
    int32_t ix = 0;

    // The first line tells us whether the file is sorted

    std::string line = block->str();
    header.sorted    = _is_sorted(line);
    ix += static_cast<int32_t>(line.length() + 1);

    // The next set of lines are for the reference sequences

    while (ix < length) {

        line = block->str();
        ix += static_cast<int32_t>(line.length() + 1);

        // If we have reached the aux data, we are done with the references

        if (_is_aux(line)) { break; }

        // Else, increment the number of references

        if (_is_seq(line)) { header.references++; }

    }

    // Skip any remaining blank space

    if (ix < length) { (void)block->bytes(length - ix); }

    return header;

}


Header _read_sam_header(std::unique_ptr<ByteStream>& block) {

    return _read_sam_header(block, block->remaining());

}





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


CIGAR_op::CIGAR_op(CIGAR_t type, int32_t length) : _type(type), _length(length) {}


CIGAR_t CIGAR_op::type() const {

    return _type;

}


base_t CIGAR_op::last() const {

    return base;

}


base_t CIGAR_op::last(base_t reference) const {

    if (!sub.empty()) {
        return sub[reference];
    } else {
        return base;
    }

}


int32_t CIGAR_op::length() const {

    return _length;

}


int32_t CIGAR_op::rlength() const {

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


int32_t CIGAR_op::qlength() const {

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


void CIGAR_op::extend(int32_t val) {

    _length += val;

}


bool CIGAR_op::empty() const {

    return _length == 0;

}


bool CIGAR_op::advance(int32_t offset) {

    if (_pos < _length + offset) {
        _pos++;
        return true;
    }
    return false;

}


std::pair<CIGAR_op, CIGAR_op> CIGAR_op::split(int32_t n) const {

    int32_t l1 = std::min(_length, n);
    int32_t l2 = _length - l1;

    return {CIGAR_op(_type, l1), CIGAR_op(_type, l2)};

}


CIGAR_op CIGAR_op::match() const {

    return {CIGAR_t::MATCH, _length};

}


CIGAR::CIGAR(int32_t size) {

    _str.reserve(size);

}


bool CIGAR::empty() const {

    return _str.empty();

}


int32_t CIGAR::rlength() const {

    int32_t _rlength = 0;
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


void CIGAR::extend(const CIGAR_op& op) {

    if (size() > 0 && _str.back().type() == op.type()) {
        _str.back().extend(op.length());
        _str.back().base = op.base;
        _str.back().sub  = op.sub;
    } else {
        _str.push_back(op);
    }

}


void CIGAR::append(const CIGAR_op& op) {

    _str.push_back(op);

}


int32_t CIGAR::size() const {

    return static_cast<int32_t>(_str.size());

}


int32_t CIGAR::bases() const {

    int32_t bases = 0;
    for (const auto& op : _str) {
        bases += op.length();
    }

    return bases;

}


CIGAR_op CIGAR::back() const {

    if (!_str.empty()) { return _str.back(); }
    return {};

}


auto CIGAR::begin() -> decltype(_str.begin()) { return _str.begin(); }
auto CIGAR::begin() const -> decltype(_str.begin()) { return _str.begin(); }
auto CIGAR::end() -> decltype(_str.end()) { return _str.end(); }
auto CIGAR::end() const -> decltype(_str.end()) { return _str.end(); }





//
// PHRED
//





template <typename dtype>
static inline std::vector<dtype> _get_mask(
    std::span<const qual_t> quality,
    qual_t min,
    int32_t window
) {

    auto length = static_cast<int32_t>(quality.size());
    int32_t bad = 0;
    window = std::min(window, length);

    int32_t start = 0;
    int32_t stop  = 0;

    if (length < 2 * window + 1) {
        start = length - window;
        stop  = window + 1;
    } else {
        start = window + 1;
        stop  = length - window;
    }

    std::vector<dtype> mask(length + 1, 0);

    // 1. Count the bad bases in the starting window region
    for (int32_t ix = 0; ix < window; ix++) {
        bad += (quality[ix] < min);
    }

    // 2. Increment bad bases at the 3' end as the window region grows
    for (int32_t ix = 0; ix < start; ix++) {
        bad += (quality[ix + window] < min);
        mask[ix] = (bad == 0);
    }

    // 3. Increment and decrement bad bases as the window region shifts
    if (length < 2 * window + 1) {
        for (int32_t ix = start; ix < stop; ix++) {
            mask[ix] = (bad == 0);
        }
    } else {
        for (int32_t ix = start; ix < stop; ix++) {
            bad -= (quality[ix - window - 1] < min);
            bad += (quality[ix + window] < min);
            mask[ix] = (bad == 0);
        }
    }

    // 4. Decrement bad bases at the 5' end as the window region shrinks
    for (int32_t ix = stop; ix < length + 1; ix++) {
        bad -= (quality[ix - window - 1] < min);
        mask[ix] = (bad == 0);
    }

    return mask;

}


static inline std::string _phred_to_str(const std::vector<qual_t>& qualities) {

    auto _length = static_cast<int32_t>(qualities.size());
    std::string str(_length, '\0');

    for (int32_t ix = 0; ix < _length; ix++) {
        str[ix] = static_cast<char>(qualities[ix] + PHRED_OFFSET);
    }

    return str;

}


PHRED::PHRED(const std::vector<qual_t>& qualities) : _qualities(qualities) {}


qual_t PHRED::operator[](int32_t ix) const {

    return _qualities[ix];

}

std::string PHRED::str() const {

    return _phred_to_str(_qualities);

}


bool PHRED::check(int32_t ix, qual_t min, int32_t window) const {

    auto _length = static_cast<int32_t>(_qualities.size());
    auto _zero   = static_cast<int32_t>(0);

    int32_t jx_min = std::max(ix - window, _zero);
    int32_t jx_max = std::min(ix + window + 1, _length);
    bool quality = true;

    for (int32_t jx = jx_min; jx < jx_max; jx++) {
        quality &= (_qualities[jx] >= min);
        quality &= (_qualities[jx] < MISSING_MAPQ);
    }

    return quality;

}


template <typename dtype>
std::vector<dtype> PHRED::mask(qual_t min, int32_t window) const {

    return _get_mask<dtype>(_qualities, min, window);

}





//
// AlignmentIter
//





Iterator::Iterator(int64_t reads) : _reads(reads) {}


bool Iterator::end() const noexcept {

    return _curr >= _reads;

}


Alignment EmptyIterator::next() {

    throw std::runtime_error("Cannot get alignments from an empty iterator.");

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


void IndexBlock::write_bad_ptr(std::ofstream& file) {

    int64_t _tmp_ptr = -1;
    file.write(reinterpret_cast<char*>(&_tmp_ptr), sizeof(int64_t));

}


void IndexBlock::write_reads(std::ofstream& file) {

    file.write(reinterpret_cast<char*>(&reads), sizeof(int64_t));

}


Index::Index(const std::string& filename, int32_t references)
    : _name(filename), _file(filename), _references(references) {}


IndexBlock Index::read(int32_t ix) {

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
    for (int32_t ix = 0; ix < _references; ix++) {
        reads += read(ix).reads;
    }

    return reads;

}


int64_t Index::unaligned() {

    int64_t unaligned = 0;

    _file.seekg(_references * sizeof(IndexBlock));
    _file.read(reinterpret_cast<char*>(&unaligned), sizeof(int64_t));

    return unaligned;

}





//
// File
//





File::File(const std::string& name, FileType type)
    : _name(name), _type(type), _hts_bgzf(_open_bgzf(name)) {}


File::~File() {

    _close_bgzf(_hts_bgzf);

}


std::string File::name() const {

    return _name;

}


FileType File::type() const {

    return _type;
}


int32_t File::size() const {

    return _references;

}


int64_t File::aligned() {

    return _index.aligned();

}


int64_t File::unaligned() {

    return _index.unaligned();

}


static inline std::unique_ptr<File> _get_file(const std::string& name) {

    BGZF* file    = _open_bgzf(name);
    FileType type = _get_filetype(file);
    _close_bgzf(file);

    switch (type) {

        case FileType::SAM: {
            return _get_sam(name);
        }

        case FileType::BAM: {
            return _get_bam(name);
        }

        case FileType::CRAM: {
            return _get_cram(name);
        }

    }

}





//
// FileGroup
//





FileGroup::FileGroup(const std::vector<std::string>& filenames) {

    _throw_if_has_duplicate_paths(filenames);

    for (const auto& name : filenames) {
        _group.push_back(_get_file(name));
    }

}


int32_t FileGroup::size() const {

    return static_cast<int32_t>(_group.size());

}


int64_t FileGroup::aligned() {

    int64_t _reads = 0;
    for (auto& file : _group) {
        _reads += file->aligned();
    }

    return _reads;

}


int64_t FileGroup::unaligned() {

    int64_t _reads = 0;
    for (auto& file : _group) {
        _reads += file->unaligned();
    }

    return _reads;

}


int32_t FileGroup::references() const {

    if (size() > 0) {
        return _group[0]->size();
    } else {
        return 0;
    }

}


auto FileGroup::begin() -> decltype(_group.begin()) { return _group.begin(); }


auto FileGroup::begin() const -> decltype(_group.begin()) { return _group.begin(); }


auto FileGroup::end() -> decltype(_group.end()) { return _group.end(); }


auto FileGroup::end() const -> decltype(_group.end()) { return _group.end(); }





//
// Printing utilities
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


std::ostream& operator<<(std::ostream& os, const CIGAR_op& op) {

    os << std::left << std::setw(9) << op.type() << std::right << std::setw(4) << op.length();

    return os;

}



std::ostream& operator<<(std::ostream& os, const CIGAR& cigar) {

    for (const auto& op : cigar) {
        os << op  << "\n";
    }

    return os;

}





//
// Templates
//




template std::vector<float> PHRED::mask(qual_t min, int32_t window) const;





} // namespace HTS


