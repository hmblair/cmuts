#include "fasta.hpp"





//
// Tools for working with binary sequence data
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


static inline int32_t _bytes_from_seq(int32_t n) {

    return (n + BASES_PER_BYTE - 1) / BASES_PER_BYTE;

}


static inline base_t _binary_mask(base_t nuc, int32_t i) {

    return static_cast<base_t>(
        nuc << ((3 - (i % BASES_PER_BYTE)) * 2)
    );

}


static inline seq_t _to_binary(const std::string& sequence) {

    int32_t len   = sequence.length();
    int32_t bytes = _bytes_from_seq(len);
    seq_t binary(bytes);

    for (int32_t i = 0; i < len; i++) {

        base_t nuc = _to_int(sequence[i]);
        if (nuc == IX_UNK) {
            throw std::runtime_error("Invalid nucleotide \"" + std::to_string(sequence[i]) + "\" encountered.");
        }

        binary[i / BASES_PER_BYTE] |= _binary_mask(nuc, i);

    }

    return binary;

}


static inline base_t _from_binary(base_t binary, int32_t j) {

    return (binary >> ((3 - j) * 2)) & BIN_T;

}


static inline int32_t _mod_round_up(int32_t val, int32_t div) {

    int32_t out = val % div;
    if (out == 0) { out = div; }
    return out;

}


static inline seq_t _from_binary(const seq_t& binary, int32_t length) {

    int32_t bytes = binary.size();
    seq_t sequence(length);

    int32_t ix, jx;
    for (ix = 0; ix < bytes - 1; ix++) {
        for (jx = 0; jx < BASES_PER_BYTE; jx++) {
            sequence[(BASES_PER_BYTE * ix) + jx] = _from_binary(binary[ix], jx);
        }
    }

    int32_t end = _mod_round_up(length, BASES_PER_BYTE);
    for (jx = 0; jx < end; jx++) {
        sequence[(BASES_PER_BYTE * ix) + jx] = _from_binary(binary[ix], jx);
    }

    return sequence;

}


static inline seq_t _read_binary_sequence(std::ifstream& file, int32_t length, int64_t offset) {

    int32_t bytes = _bytes_from_seq(length);
    seq_t binary(bytes);

    file.seekg(offset);
    file.read(reinterpret_cast<char*>(binary.data()), bytes * sizeof(base_t));

    return _from_binary(binary, length);

}


static inline void _write_binary_sequence(std::ofstream& file, const std::string& sequence) {

    seq_t binary = _to_binary(sequence);
    file.write(reinterpret_cast<const char*>(binary.data()), binary.size() * sizeof(base_t));

}


static inline HeaderBlock _read_fasta_header_block(std::ifstream& file) {

    HeaderBlock block;
    file.read(reinterpret_cast<char*>(&block.length), sizeof(int32_t));
    file.read(reinterpret_cast<char*>(&block.sequences), sizeof(int32_t));
    return block;

}


static inline std::vector<HeaderBlock> _read_fasta_header(std::ifstream& file) {

    HeaderBlock block;
    std::vector<HeaderBlock> blocks;

    while (block.length != EOH) {

        block = _read_fasta_header_block(file);
        if (!block.empty()) { blocks.push_back(block); }

    }

    return blocks;

}


static inline void _write_fasta_header_block(std::ofstream& file, const HeaderBlock& block) {

    file.write(reinterpret_cast<const char*>(&block.length), sizeof(int32_t));
    file.write(reinterpret_cast<const char*>(&block.sequences), sizeof(int32_t));

}


static inline void _write_fasta_eoh(std::ofstream& file) {

    file.write(reinterpret_cast<const char*>(&EOH), sizeof(int32_t));

}


static inline void _write_fasta_header(std::ofstream& file, const std::vector<HeaderBlock>& blocks) {

    for (const auto& block : blocks) {
        if (!block.empty()) { _write_fasta_header_block(file, block); }
    }

    _write_fasta_eoh(file);

}


static inline std::vector<HeaderBlock> _get_fasta_header(const std::string& name) {

    FASTA fasta(name);
    HeaderBlock block;
    std::vector<HeaderBlock> blocks;

    std::string sequence = fasta.next();
    while (!sequence.empty()) {

        auto length = static_cast<int32_t>(sequence.length());
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


static inline void _fasta_to_binary(const std::string& fname, const std::string& bname) {

    __log(
        _LOG_FILE,
        "Building index for " + fname + "..."
    );

    _throw_if_not_exists(fname);
    _throw_if_exists(bname);

    std::ofstream binary(bname);
    std::vector<HeaderBlock> blocks = _get_fasta_header(fname);
    _write_fasta_header(binary, blocks);

    FASTA fasta(fname);

    std::string sequence = fasta.next();
    if (sequence.empty()) {
        binary.close();
        return;
    }

    int32_t count = 0;

    _write_binary_sequence(binary, sequence);
    while (!sequence.empty()) {
        sequence = fasta.next();
        _write_binary_sequence(binary, sequence);
        count++;
    }

    binary.close();

    __log(
        _LOG_FILE,\
        "Successfully created " + bname + " with " + std::to_string(count) + " sequences."
    );

}





//
// FASTA Header
//





bool HeaderBlock::empty() const {

    return (length <= 0 || sequences == 0);

}


Header::Header(std::ifstream& file) :
    _blocks(_read_fasta_header(file)),
    _offset((2 * _blocks.size() + 1) * sizeof(int32_t)) {}


int32_t Header::size() const {

    int32_t _sequences = 0;
    for (const auto& block : _blocks) {
        _sequences += block.sequences;
    }

    return _sequences;

}


int32_t Header::length(int32_t ix) const {

    int32_t count = 0;
    for (const auto& block : _blocks) {

        if (count + block.sequences > ix) {
            return block.length;
        }
        count += block.sequences;

    }

    throw std::runtime_error("The index " + std::to_string(ix) + " is larger than the number of seqeunces in the header.");

}


int32_t Header::longest() const {

    int32_t _longest = 0;
    for (const auto& block : _blocks) {
        _longest = std::max(_longest, block.length);
    }

    return _longest;

}


Offset Header::offset(int32_t ix) {

    Offset offset;
    int32_t count  = 0;

    // Loop until we have consumed ix sequences

    for (const auto& block : _blocks) {

        auto bytes = static_cast<int64_t>(
            _bytes_from_seq(block.length)
        );

        if (count + block.sequences > ix) {
            offset.offset += (ix - count) * bytes;
            offset.length = block.length;
            break;
        } else {
            offset.offset += block.sequences * bytes;
            count += block.sequences;
        }

    }

    // Account for the header by adding _offset

    offset.offset += _offset;
    return offset;

}




//
// Interface with HDF5
//




template <typename dtype>
static inline HDF5::Memspace<dtype, FASTA_DATASET_SIZE> __memspace(
    const BinaryFASTA& fasta,
    HDF5::File& hdf5
) {
    std::vector<size_t> dims = {
        static_cast<size_t>(fasta.size()),
        static_cast<size_t>(fasta.longest())
    };
    return hdf5.memspace<dtype, FASTA_DATASET_SIZE>(dims, FASTA_DATASET);
}

template <typename dtype>
static inline void _fasta_to_hdf5(
    BinaryFASTA& fasta,
    HDF5::File& hdf5,
    const MPI::Manager& mpi
) {

    HDF5::Memspace<dtype, FASTA_DATASET_SIZE> memspace = __memspace<dtype>(fasta, hdf5);
    MPI::Chunk chunk = mpi.chunk(
        hdf5.chunk_size(),
        fasta.size()
    );

    for (int32_t ix = chunk.low; ix < chunk.high; ix+= chunk.step) {
        for (int32_t jx = 0; jx < chunk.size && ix + jx < fasta.size(); jx++) {

            view_t<dtype, FASTA_DATASET_SIZE> arr = memspace.view(jx);
            seq_t sequence = fasta.sequence(ix + jx);
            std::copy(sequence.begin(), sequence.end(), arr.begin());

        }

        memspace.safe_write(ix);
        memspace.clear();

    }

}





//
// BinaryFASTA
//



static inline std::string _cmfa_extension(const std::string& filename) {

    size_t last_dot = filename.find_last_of('.');
    std::string base = (last_dot != std::string::npos) ? filename.substr(0, last_dot) : filename;

    return base + ".cmfa";

}


BinaryFASTA::BinaryFASTA(const std::string& fasta, bool rebuild)
    : _fasta_name(fasta), _name(_cmfa_extension(fasta)) {

    _throw_if_not_exists(_fasta_name);

    if (rebuild) { _delete(_name); }

    if (!_exists(_name) && !cmuts::mutex::check(_name)) {
        cmuts::mutex::Mutex mutex = cmuts::mutex::lock(_name);
        _fasta_to_binary(_fasta_name, _name);
    }

    // Will trigger if another thread started creating the index file first

    cmuts::mutex::wait(_name);

    _file   = std::ifstream(_name);
    _header = Header(_file);

    __log(_LOG_FILE, "Successfully loaded " + _name + ".");

}


BinaryFASTA::BinaryFASTA(BinaryFASTA&& other) noexcept
    : _name(std::move(other._name)),
      _file(std::move(other._file)),
      _header(std::move(other._header)) {}


BinaryFASTA& BinaryFASTA::operator=(BinaryFASTA&& other) noexcept {

    if (this != &other) {
        _name = other._name;
        _file = std::move(other._file);
        _header = std::move(other._header);
    }

    return *this;

}


std::string BinaryFASTA::name() const {

    return _name;

}


int32_t BinaryFASTA::size() const {

    return _header.size();

}


int32_t BinaryFASTA::length(int32_t ix) const {

    return _header.length(ix);

}


seq_t BinaryFASTA::sequence(int32_t ix) {

    Offset offset = _header.offset(ix);

    return _read_binary_sequence(_file, offset.length, offset.offset);

}


int32_t BinaryFASTA::longest() const {

    return _header.longest();

}


void BinaryFASTA::hdf5(HDF5::File& hdf5, const MPI::Manager& mpi) {

    _fasta_to_hdf5<base_t>(*this, hdf5, mpi);

}
