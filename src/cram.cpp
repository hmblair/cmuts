#include "cram.hpp"





namespace HTS {



static inline base_t _ascii_to_int(uint8_t base) {

    switch (base) {
        case static_cast<uint8_t>(A):   { return IX_A;   }
        case static_cast<uint8_t>(C):   { return IX_C;   }
        case static_cast<uint8_t>(G):   { return IX_G;   }
        case static_cast<uint8_t>(T):   { return IX_T;   }
        case static_cast<uint8_t>(U):   { return IX_U;   }
        default:                        { return IX_UNK; }
    }

}


template <typename dtype>
static inline int32_t _find_first_occurence(
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



static inline int32_t _read_ITF8(BGZF* file) {

    // Read the first byte, which will tell us the size of the ITF8

    static uint8_t buffer[5];
    _read_bgzf(file, buffer, 1);

    // Compute the size of the ITF8 by counting the leading ones, up to a max of 4

    int32_t size = std::countl_one(static_cast<uint8_t>(buffer[0] & 0xF0));

    // Initialize the value with the least 7 - size significant bits

    int32_t value = static_cast<int32_t>(buffer[0] & ITF8_MASK[size]);
    if (size == 0) { return value; }

    // Read the appropriate amount of bytes

    _read_bgzf(file, buffer + 1, size);

    // Compute the value by appending the data from the read bytes

    for (int32_t ix = 1; ix < size; ix++) {
        value = (value << 8) | buffer[ix];
    }

    if (size == 4) {
        value = (value << 4) | (buffer[4] & 0x0f);
    } else {
        value = (value << 8) | buffer[size];
    }

    return value;

}


static inline int64_t _read_LTF8(BGZF* file) {

    static uint8_t buffer[9];
    _read_bgzf(file, buffer, 1);

    int32_t size = std::countl_one(buffer[0]);
    if (size > 0) { _read_bgzf(file, buffer + 1, size); }

    int64_t value = static_cast<int64_t>(buffer[0]);
    for (int32_t ix = 1; ix <= size; ix++) {
        value = (value << 8) | buffer[ix];
    }

    return value;

}


static inline std::vector<int32_t> _read_ITF8_array(BGZF* _bgzf_file) {

    int32_t size = _read_ITF8(_bgzf_file);
    std::vector<int32_t> array(size);

    for (int32_t ix = 0; ix < size; ix++) {
        array[ix] = _read_ITF8(_bgzf_file);
    }

    return array;

}


static inline ITF8 _decode_ITF8(const uint8_t* buffer, int32_t remaining) {

    ITF8 itf8;

    // Compute the size of the ITF8 by counting the leading ones, up to a max of 4

    int32_t bytes = std::countl_one(static_cast<uint8_t>(buffer[0] & 0xF0));
    itf8.size += bytes;

    // Initialize the value with the least 7 - size significant bits

    itf8.value = static_cast<int32_t>(buffer[0] & ITF8_MASK[bytes]);

    // Only read the bytes that are available in the buffer

    bytes = std::min(bytes, remaining - 1);
    itf8.read += bytes;

    if (bytes == 0) { return itf8; }

    // Compute the value by appending the data from the read bytes

    for (int32_t ix = 1; ix < bytes; ix++) {
        itf8.value = (itf8.value << 8) | buffer[ix];
    }

    if (bytes == 4) {
        itf8.value = (itf8.value << 4) | (buffer[bytes] & 0x0F);
    } else {
        itf8.value = (itf8.value << 8) | buffer[bytes];
    }

    return itf8;

}


static inline LTF8 _decode_LTF8(const uint8_t* buffer, int32_t remaining) {

    LTF8 ltf8;

    // Similar to _decode_ITF8, just with 8 size bits rather than 4

    int32_t bytes = std::countl_one(static_cast<uint8_t>(buffer[0]));
    ltf8.size += bytes;
    ltf8.value = static_cast<int64_t>(buffer[0]);

    // Only read the bytes that are available in the buffer

    bytes = std::min(bytes, remaining - 1);
    ltf8.read += bytes;

    if (bytes == 0) { return ltf8; }

    for (int32_t ix = 1; ix <= bytes; ix++) {
        ltf8.value = (ltf8.value << 8) | buffer[ix];
    }

    return ltf8;

}


static inline std::unique_ptr<ByteStream> _get_stream(const cramBlock& block) {

    std::span data(block.data.begin(), block.data.end());

    switch (block.method) {

        case CompressionMethod::None:    {
            return std::make_unique<dataStream>(data);
        }

        case CompressionMethod::gzip:    {
            return std::make_unique<zlibStream>(data, block.raw, BGZF_BUFFER);
        }

        case CompressionMethod::rans4x8: {
            return std::make_unique<ransStream>(data, block.raw);
        }

        default: {
            throw std::runtime_error("Unsupported compression method \"" + std::to_string(static_cast<int32_t>(block.method)) + "\".");
        }

    }

}


static inline std::unique_ptr<ByteStream> _get_file_stream(const cramBlock& block) {

    switch (block.method) {

        case CompressionMethod::None:    {
            return std::make_unique<bgzfFileStream>(block.file, block.size);
        }

        case CompressionMethod::gzip:    {
            return std::make_unique<zlibFileStream>(block.file, block.size, block.raw, BGZF_BUFFER);
        }

        default: {
            throw std::runtime_error("Unsupported file-streaming compression method \"" + std::to_string(static_cast<int32_t>(block.method)) + "\".");
        }

    }

}










template <typename dtype>
static inline std::vector<dtype> _read_bgzf_array(BGZF* _bgzf_file) {

    int32_t count = _read_ITF8(_bgzf_file);

    std::vector<dtype> values(count);
    _read_bgzf(_bgzf_file, values.data(), sizeof(dtype) * count);

    return values;

}


static inline Codec_t _codec_from_cram(int32_t value) {

    switch (value) {

        case 1:  { return Codec_t::ID;     }
        case 3:  { return Codec_t::Huffman;  }
        case 4:  { return Codec_t::ByteArrayLen; }
        case 5:  { return Codec_t::ByteStop; }
        case 6:  { return Codec_t::Beta;     }
        case 7:  { return Codec_t::SubExp;     }
        default: {
            throw std::runtime_error("Invalid codec \"" + std::to_string(value) + "\".");
        }

    }

}



static inline CompressionMethod _get_compression(uint8_t value) {

    switch (value) {

        case 0:  { return CompressionMethod::None;    }
        case 1:  { return CompressionMethod::gzip;    }
        case 4:  { return CompressionMethod::rans4x8; }
        default: {
            throw std::runtime_error("Invalid compression method \"" + std::to_string(value) + "\".");
        }

    }

}


int64_t cramContainer::aligned() {

    int64_t reads = 0;
    for (auto& slice: slices) {
        reads += slice.aligned();
    }
    return reads;

}


int64_t cramContainer::unaligned() {

    int64_t reads = 0;
    for (auto& slice: slices) {
        reads += slice.unaligned();
    }
    return reads;

}


static inline ExtData_t _aux_from_code(const std::string &code) {

    if (code == "BF") { return ExtData_t::BF; }
    if (code == "CF") { return ExtData_t::CF; }
    if (code == "RI") { return ExtData_t::RI; }
    if (code == "RL") { return ExtData_t::RL; }
    if (code == "AP") { return ExtData_t::AP; }
    if (code == "RG") { return ExtData_t::RG; }
    if (code == "RN") { return ExtData_t::RN; }
    if (code == "MF") { return ExtData_t::MF; }
    if (code == "NS") { return ExtData_t::NS; }
    if (code == "NP") { return ExtData_t::NP; }
    if (code == "TS") { return ExtData_t::TS; }
    if (code == "NF") { return ExtData_t::NF; }
    if (code == "TL") { return ExtData_t::TL; }
    if (code == "FN") { return ExtData_t::FN; }
    if (code == "FC") { return ExtData_t::FC; }
    if (code == "FP") { return ExtData_t::FP; }
    if (code == "DL") { return ExtData_t::DL; }
    if (code == "BB") { return ExtData_t::BB; }
    if (code == "QQ") { return ExtData_t::QQ; }
    if (code == "BS") { return ExtData_t::BS; }
    if (code == "IN") { return ExtData_t::IN; }
    if (code == "RS") { return ExtData_t::RS; }
    if (code == "PD") { return ExtData_t::PD; }
    if (code == "HC") { return ExtData_t::HC; }
    if (code == "SC") { return ExtData_t::SC; }
    if (code == "MQ") { return ExtData_t::MQ; }
    if (code == "BA") { return ExtData_t::BA; }
    if (code == "QS") { return ExtData_t::QS; }

    throw std::invalid_argument("Unknown auxilliary CRAM code " + code + ".");

}

static inline std::string _code_from_aux(ExtData_t aux) {

    if (aux == ExtData_t::BF) { return "BF"; }
    if (aux == ExtData_t::CF) { return "CF"; }
    if (aux == ExtData_t::RI) { return "RI"; }
    if (aux == ExtData_t::RL) { return "RL"; }
    if (aux == ExtData_t::AP) { return "AP"; }
    if (aux == ExtData_t::RG) { return "RG"; }
    if (aux == ExtData_t::RN) { return "RN"; }
    if (aux == ExtData_t::MF) { return "MF"; }
    if (aux == ExtData_t::NS) { return "NS"; }
    if (aux == ExtData_t::NP) { return "NP"; }
    if (aux == ExtData_t::TS) { return "TS"; }
    if (aux == ExtData_t::NF) { return "NF"; }
    if (aux == ExtData_t::TL) { return "TL"; }
    if (aux == ExtData_t::FN) { return "FN"; }
    if (aux == ExtData_t::FC) { return "FC"; }
    if (aux == ExtData_t::FP) { return "FP"; }
    if (aux == ExtData_t::DL) { return "DL"; }
    if (aux == ExtData_t::BB) { return "BB"; }
    if (aux == ExtData_t::QQ) { return "QQ"; }
    if (aux == ExtData_t::BS) { return "BS"; }
    if (aux == ExtData_t::IN) { return "IN"; }
    if (aux == ExtData_t::RS) { return "RS"; }
    if (aux == ExtData_t::PD) { return "PD"; }
    if (aux == ExtData_t::HC) { return "HC"; }
    if (aux == ExtData_t::SC) { return "SC"; }
    if (aux == ExtData_t::MQ) { return "MQ"; }
    if (aux == ExtData_t::BA) { return "BA"; }
    if (aux == ExtData_t::QS) { return "QS"; }

    return "NULL";

}


static inline AuxDataEncoding _get_aux_enc_type(ExtData_t type) {

    switch (type) {

        case ExtData_t::BF:
        case ExtData_t::CF:
        case ExtData_t::RI:
        case ExtData_t::RL:
        case ExtData_t::AP:
        case ExtData_t::RG:
        case ExtData_t::MF:
        case ExtData_t::NS:
        case ExtData_t::NP:
        case ExtData_t::TS:
        case ExtData_t::NF:
        case ExtData_t::TL:
        case ExtData_t::FN:
        case ExtData_t::FP:
        case ExtData_t::DL:
        case ExtData_t::RS:
        case ExtData_t::PD:
        case ExtData_t::HC:
        case ExtData_t::MQ: {

            return AuxDataEncoding::Int32;

        }

        case ExtData_t::FC:
        case ExtData_t::BS:
        case ExtData_t::BA:
        case ExtData_t::QS: {

            return AuxDataEncoding::Byte;

        }

        case ExtData_t::RN:
        case ExtData_t::BB:
        case ExtData_t::QQ:
        case ExtData_t::IN:
        case ExtData_t::SC: {

            return AuxDataEncoding::ByteArray;

        }

    }

}


Codec::Codec(Codec_t type)
    : _type(type) {}


Codec::Codec(Codec_t type, int32_t id)
    : _type(type), _id(id) {}


ExternalCodec::ExternalCodec(int32_t id)
    : Codec(Codec_t::ID, id) {}


ByteArrayStopCodec::ByteArrayStopCodec(int32_t id, uint8_t stop)
    : Codec(Codec_t::ByteStop, id), _stop(stop) {}


HuffmanCodec::HuffmanCodec(const std::vector<int32_t>& alphabet, const std::vector<int32_t>& lengths)
    : Codec(Codec_t::Huffman), _alphabet(alphabet), _lengths(lengths) {}


BetaCodec::BetaCodec(int32_t offset, int32_t length)
    : Codec(Codec_t::Beta), _offset(offset), _length(length) {}


SubExpCodec::SubExpCodec(int32_t offset, int32_t order)
    : Codec(Codec_t::SubExp), _offset(offset), _order(order) {}


Codec_t Codec::type() const { return _type; }

int32_t Codec::id() const { return _id; }

void Codec::attach(const std::shared_ptr<ByteStream>& stream) { _stream = stream; }


// Default implementations


uint8_t Codec::byte() { return _stream->byte(); }


int32_t Codec::integer() { return _stream->itf8(); }


std::vector<uint8_t> Codec::array() { 

    throw std::runtime_error("No default for array().");

}

std::vector<uint8_t> Codec::array(int32_t length) { 

    throw std::runtime_error("No default for array(int32_t).");

}


// Codec-specific implementations


std::vector<uint8_t> ExternalCodec::array(int32_t length) { return _stream->bytes(length); }


std::vector<uint8_t> ByteArrayStopCodec::array() { return _stream->line(_stop); }


uint8_t HuffmanCodec::byte() {

    return static_cast<uint8_t>(_alphabet[0]);

}


int32_t HuffmanCodec::integer() {

    return _alphabet[0];

}


int32_t BetaCodec::integer() {

    return _stream->bits(_length) + _offset;

}


int32_t SubExpCodec::integer() {

    // Read the number of leading ones.
    // This also takes care of the trailing zero.

    uint8_t bit   = _stream->bits(1);
    uint8_t count = 0;

    while (bit == 1) {
        count++;
        bit = _stream->bits(1);
    }

    // Decode as per the standard

    int32_t length = 0;

    if (count == 0) {
        length = static_cast<int32_t>(_order);
    } else {
        length = static_cast<int32_t>(count + _order - 1);
    }

    int32_t value = static_cast<int32_t>(_stream->bits(length));
    if (count > 0) { value += (1L << length); }

    return value - _offset;

}





//
// cramBlock
//





cramBlock::cramBlock(BGZF* file)
    : method  (_get_compression(_read_bgzf_single<uint8_t>(file))),
      type    (static_cast<cramBlockType>(_read_bgzf_single<uint8_t>(file))),
      content (_read_ITF8(file)),
      size    (_read_ITF8(file)),
      raw     (_read_ITF8(file)),
      file    (file) {}


cramBlock::cramBlock(uint8_t*& _data) {

    // The compression method and block type

    method = _get_compression(*_data);           _data++;
    type   = static_cast<cramBlockType>(*_data); _data++;

    // The external content identifier, compressed data size, and
    // uncompressed data size

    ITF8 _content = _decode_ITF8(_data, ITF8_MAX_SIZE); _data += _content.read;
    ITF8 _size    = _decode_ITF8(_data, ITF8_MAX_SIZE); _data += _size.read;
    ITF8 _raw     = _decode_ITF8(_data, ITF8_MAX_SIZE); _data += _raw.read;

    content = _content.value;
    size    = _size.value;
    raw     = _raw.value;

    // The pointer to the start of the block data. Skip the CRC32 checksum
    // as well, which we have no use for

    data = std::span<const uint8_t>(_data, _data + size);
    _data += (size + CRC_SIZE);

}


void cramBlock::memcpy(uint8_t* dest, int32_t size) {

    std::memcpy(dest, data.data(), size);
    data = data.subspan(size);

}


uint8_t cramBlock::byte() {

    uint8_t out = data[0];
    data = data.subspan(1);
    return out;

}


std::vector<uint8_t> cramBlock::bytes(int32_t n) {

    std::vector<uint8_t> bytes(data.begin(), data.begin() + n);
    data = data.subspan(n);
    return bytes;

}


template <typename dtype>
std::vector<dtype> cramBlock::array() {

    int32_t size = itf8();
    std::vector<uint8_t> _bytes = bytes(size * sizeof(dtype));

    std::vector<dtype> result(size);
    std::memcpy(result.data(), _bytes.data(), size * sizeof(dtype));

    return result;

}


int32_t cramBlock::itf8() { 

    ITF8 _itf8 = _decode_ITF8(data.data(), data.size());
    data = data.subspan(_itf8.read);
    return _itf8.value;

}


int64_t cramBlock::ltf8() {

    LTF8 _ltf8 = _decode_LTF8(data.data(), data.size());
    data = data.subspan(_ltf8.read);
    return _ltf8.value;

}


std::vector<int32_t> cramBlock::itf8_arr() {

    // The first ITF8 value is the array size

    int32_t n = itf8();
    std::vector<int32_t> values(n);

    // The remaining ITF8 values are the array elements

    for (int32_t ix = 0; ix < n; ix++) { values[ix] = itf8(); }

    return values;

}


std::vector<int64_t> cramBlock::ltf8_arr() {

    // The first ITF8 value is the array size

    int32_t n = itf8();
    std::vector<int64_t> values(n);

    // The remaining LTF8 values are the array elements

    for (int32_t ix = 0; ix < n; ix++) { values[ix] = ltf8(); }

    return values;

}


std::string cramBlock::str(int32_t length) {

    std::string result(data.begin(), data.begin() + length);
    data = data.subspan(length);
    return result;

}


static inline std::vector<std::vector<base_t>> _parse_substitution_matrix(uint8_t* array) {

    std::vector<std::vector<base_t>> decoder(BASES + 1);

    for (base_t ix = 0; ix < BASES; ix++) {

        uint8_t line   = array[ix];
        uint8_t offset = 0;

        for (base_t jx = 0; jx < BASES + 1; jx++) {

            if (ix == jx) { continue; }

            if (jx >= ix) {
                offset = 1;
            } else {
                offset = 0;
            }

            uint8_t code = line >> (6 - 2 * (jx - offset));
            code &= 0x3;

            decoder[code].push_back(jx); 


        }

    }

    return decoder;

}

SubstitutionMatrix::SubstitutionMatrix(uint8_t* matrix)
    : _decoder(_parse_substitution_matrix(matrix)) {}


std::vector<base_t> SubstitutionMatrix::parse(base_t code) { return _decoder[code]; }


static inline PreservationMap _get_preservation_map(cramBlock& block) {

    PreservationMap map;

    int32_t size = block.itf8();
    int32_t keys = block.itf8();

    for (int32_t ix = 0; ix < keys; ix++) {

        std::string key = block.str(KEY_SIZE);

        if (key == "TD") {
            map.TD = block.array<uint8_t>();
        } else if (key == "SM") {
            block.memcpy(map.SM, SM_SIZE);
        } else if (key == "RR") {
            map.RR = static_cast<bool>(block.byte());
        } else if (key == "AP") {
            map.AP = static_cast<bool>(block.byte());
        } else if (key == "RN") {
            map.RN = static_cast<bool>(block.byte());
        } else {
            throw std::runtime_error("Invalid preservation map key \"" + key + "\".");
        }

    }

    return map;

}


static inline ExtData_t _fc_to_enum(uint8_t value) {

    switch (value) {

        case 0x62: { return ExtData_t::BB; }
        case 0x71: { return ExtData_t::QQ; }
        case 0x42: { return ExtData_t::QS; }
        case 0x58: { return ExtData_t::BS; }
        case 0x49: { return ExtData_t::IN; }
        case 0x44: { return ExtData_t::DL; }
        case 0x69: { return ExtData_t::BA; }
        case 0x51: { return ExtData_t::QS; }
        case 0x4E: { return ExtData_t::RS; }
        case 0x53: { return ExtData_t::SC; }
        case 0x50: { return ExtData_t::PD; }
        case 0x48: { return ExtData_t::HC; }
        default:   {
            throw std::runtime_error("Invalid value \"" + std::to_string(value) + "\" in _fc_to_enum.");
        }

    }

}



static inline std::shared_ptr<Codec> _get_codec(cramBlock& block) {

    Codec_t codec = _codec_from_cram(block.itf8());
    int32_t size  = block.itf8();

    switch (codec) {

        case Codec_t::ID: {

            int32_t id = block.itf8();
            return std::make_shared<ExternalCodec>(id);

        }

        case Codec_t::Huffman: {

            std::vector<int32_t> alphabet = block.itf8_arr();
            std::vector<int32_t> lengths  = block.itf8_arr();
            return std::make_shared<HuffmanCodec>(alphabet, lengths);

        }

        case Codec_t::ByteArrayLen: {

            (void)_get_codec(block);
            return _get_codec(block);

        }

        case Codec_t::ByteStop: {

            uint8_t stop = block.byte();
            int32_t id   = block.itf8();
            return std::make_shared<ByteArrayStopCodec>(id, stop);

        }

        case Codec_t::Beta: {

            int32_t offset = block.itf8();
            int32_t length = block.itf8();
            return std::make_shared<BetaCodec>(offset, length);

        }

        case Codec_t::SubExp: {

            int32_t offset = block.itf8();
            int32_t order  = block.itf8();
            return std::make_shared<SubExpCodec>(offset, order);

        }

    }

}


static inline CodecMap _read_data_series(cramBlock& block) {

    int32_t size = block.itf8();
    int32_t keys = block.itf8();
    CodecMap codecs;

    for (int32_t ix = 0; ix < keys; ix++) {

        std::string key  = block.str(KEY_SIZE);
        ExtData_t type   = _aux_from_code(key);
        auto codec = _get_codec(block);
        codecs.map[type] = codec;

        // Fill in the externals map if appropriate

        if (codec->type() == Codec_t::ID || codec->type() == Codec_t::ByteStop) {

            int32_t id = codec->id();
            codecs.externals[id] = type;

        }

    }

    return codecs;

}


static inline void _read_tag_dictionary(cramBlock& block) {

    int32_t size = block.itf8();
    int32_t keys = block.itf8();

    for (int32_t ix = 0; ix < keys; ix++) {

        int32_t key = block.itf8();
        ExtData_t type = ExtData_t::DL;
        (void)_get_codec(block);

    }

    // Function doesn't seem to do anything, but it is in the spec

}


cramSlice::cramSlice(uint8_t*& data) {

    // The slice header, which is contained in its own block 

    cramBlock block(data);
    if (block.type != cramBlockType::Slice) {
        auto _type = static_cast<int32_t>(block.type);
        throw std::runtime_error("Invalid slice header block type \"" + std::to_string(_type) + "\".");
    }

    reference  = block.itf8();
    start      = block.itf8();
    span       = block.itf8();
    records    = block.itf8();
    counter    = block.ltf8();
    nblocks    = block.itf8();

    // It is pointless to store the IDs array as the block already know their
    // own ID

    (void)block.itf8_arr();

    emb_ref_id = block.itf8();
    block.memcpy(md5, MD5_SIZE);

}


cramSlice::cramSlice(uint8_t*& data, CodecMap _codecs) 
    : cramSlice(data) {

    codecs = _codecs;

    for (int32_t ix = 0; ix < nblocks; ix++ ) {

        cramBlock block(data);

        // Attach core blocks to all codecs with an ID of -1

        if (block.type == cramBlockType::Core) {

            std::shared_ptr<ByteStream> stream = _get_stream(block);
            for (auto& codec : codecs.map) {

                if (codec.second->id() == -1) {
                    codec.second->attach(stream);
                }

            }

        }

        // Attach all remaining external blocks to their associated codec


        else {

            int32_t id = block.content;
            if (!codecs.externals.contains(id)) { continue; }
            ExtData_t type = codecs.externals.at(id);

            if (EXT_SKIP.contains(type)) { continue; }

            std::shared_ptr<ByteStream> stream = _get_stream(block);
            codecs.map.at(type)->attach(stream);

        }

    }

}


int32_t cramSlice::index() {

    if (multi()) {
        return codecs.map.at(ExtData_t::RI)->integer();
    } else {
        return reference;
    }

}


bool cramSlice::multi() const {

    return reference == -2;

}


int64_t cramSlice::aligned() {

    if (reference == UNALIGNED) { return 0; }

    if (multi()) {

        int64_t reads = 0;

        for (int32_t ix = 0; ix < records; ix++) {
            reads += static_cast<int64_t>(index());
        }

        return reads;

    }

    return static_cast<int64_t>(records);

}


int64_t cramSlice::unaligned() {

    if (reference >= 0) { return 0; }

    if (multi()) {

        int64_t reads = 0;

        for (int32_t ix = 0; ix < records; ix++) {
            reads += static_cast<int64_t>(index());
        }

        return reads;

    }

    return static_cast<int64_t>(records);

}


CompressionHeader::CompressionHeader(uint8_t*& data) {

    cramBlock block(data);
    if (block.type != cramBlockType::Compression) {
        throw std::runtime_error("Invalid compression header block type \"" + std::to_string(static_cast<int32_t>(block.type)) + "\".");
    }

    _map    = _get_preservation_map(block);
    _subs   = SubstitutionMatrix(static_cast<uint8_t*>(_map.SM));
    _codecs = _read_data_series(block);

    // Skip the tag dictionary

    _read_tag_dictionary(block);

}


CodecMap CompressionHeader::codecs() const { return _codecs; }


SubstitutionMatrix CompressionHeader::substitution() const { return _subs; }


bool CompressionHeader::delta() const { return _map.AP; }





//
// cramContainerBase; cramContainer; cramHeaderContainer
//





cramContainerBase::cramContainerBase(BGZF* file)
    : _file     (file),
      ptr       (bgzf_tell(file)),
      size      (_read_bgzf_single<int32_t>(file)),
      reference (_read_ITF8(file)),
      offset    (_read_ITF8(file)),
      span      (_read_ITF8(file)),
      records   (_read_ITF8(file)),
      counter   (_read_LTF8(file)),
      bases     (_read_LTF8(file)),
      blocks    (_read_ITF8(file)),
      landmarks (_read_ITF8_array(file)),
      crc32     (_read_bgzf_single<int32_t>(file)) {}


bool cramContainerBase::eof() const { 

    return (
        (offset    == CRAM_EOF_POS ) &&
        (size      == CRAM_EOF_SIZE) &&
        (reference == UNALIGNED    )
    );

}


bool cramContainerBase::multi() const {

    return reference == CRAM_MULTI;

}


cramContainer::cramContainer(BGZF* file)
    : cramContainerBase(file) {

    if (eof()) { return; }

    int32_t ix = 0;

    // The raw container data

    data.resize(size);
    _read_bgzf(file, data.data(), size);
    uint8_t* ptr = data.data();

    // The compression header, which is contained in its own block

    if (records > 0) { _comp = CompressionHeader(ptr); ix++; }

    // The slices, which contain the main blocks of data

    while (ix < blocks) {

        slices.emplace_back(ptr, _comp.codecs());
        ix += (slices.back().nblocks + 1);

    }

}


cramHeaderContainer::cramHeaderContainer(BGZF* file)
    : cramContainerBase(file) {

    cramBlock block(file);

    if (block.type != cramBlockType::Header) {
        throw std::runtime_error("Invalid file header block type \"" + std::to_string(static_cast<int32_t>(block.type)) + "\".");
    }

    stream = _get_file_stream(block);

    // Skip the first int, which stores the length of the header

    (void)stream->bytes(sizeof(int32_t));

}


void cramHeaderContainer::skip() {

    (void)_read_bgzf(_file, CRC_SIZE);

    for (int32_t ix = 1; ix < blocks; ix++) {

        cramBlock block(_file);
        (void)_read_bgzf(_file, block.size + CRC_SIZE);

    }

}


CodecMap cramContainer::codecs() const { return _comp.codecs(); }


SubstitutionMatrix cramContainer::substitution() const { return _comp.substitution(); }


bool cramContainer::delta() const { return _comp.delta(); }



static inline cramHeader _read_cram_header(BGZF* file) {

    cramHeader header;

    // The header begins with the file type and version

    FileType type = _get_filetype(file);
    if (type != FileType::CRAM) {
        throw std::runtime_error("Opening a non-CRAM file using _read_cram_header.");
    }
    auto major = _read_bgzf_single<uint8_t>(file);
    auto minor = _read_bgzf_single<uint8_t>(file);
    std::string version = std::to_string(major) + "." + std::to_string(minor);

    // Every CRAM header has an extra 20 bytes for the file ID

    (void)_read_bgzf(file, FILE_ID_SIZE);

    // The file header is in its own special container

    cramHeaderContainer container(file);

    // Read the raw SAM header

    Header _tmp_header = _read_sam_header(container.stream);

    header.sorted     = _tmp_header.sorted;
    header.references = _tmp_header.references;
    header.version    = version;

    // Read any additional padding blocks

    container.skip();

    return header;

}


static inline CIGAR_op _op_from_fc_enum(ExtData_t type, cramIterator& iter) {

    switch (type) {

        case ExtData_t::BS: {
            uint8_t code = iter.at(ExtData_t::BS)->byte();
            CIGAR_op op(CIGAR_t::MISMATCH);
            op.sub = iter.sub(code);
            return op;
        }

        case ExtData_t::BA: {
            uint8_t base = iter.at(ExtData_t::BA)->byte();
            CIGAR_op op(CIGAR_t::INS);
            op.base = _ascii_to_int(base);
            return op;
        }

        case ExtData_t::IN: {
            std::vector<uint8_t> bases = iter.at(ExtData_t::IN)->array();
            CIGAR_op op(CIGAR_t::INS, bases.size());
            op.base = _ascii_to_int(bases[bases.size() - 1]);
            return op;
        }

        case ExtData_t::DL: {
            int32_t length = iter.at(ExtData_t::DL)->integer();
            return {CIGAR_t::DEL, length};
        }

        case ExtData_t::SC: {
            std::vector<uint8_t> bases = iter.at(ExtData_t::SC)->array();
            return {CIGAR_t::SOFT, static_cast<int32_t>(bases.size())};
        }

        case ExtData_t::HC: {
            int32_t length = iter.at(ExtData_t::HC)->integer();
            return {CIGAR_t::HARD, length};
        }

        case ExtData_t::PD: {
            int32_t length = iter.at(ExtData_t::PD)->integer();
            return {CIGAR_t::PAD, length};
        }

        default: {
            throw std::runtime_error("Invalid ExtData_t \"" + _code_from_aux(type) + "\".");
        }

    }

}


cramIterator::cramIterator(BGZF* file, int32_t reads)
    : Iterator(reads), _file(file), _container(file), _remaining_in_slice(_container.slices[_slice].records) {

    int32_t _tmp = _container.slices[_slice].start;
    _offset      = std::max(_tmp - 1, 0);

}


bool cramIterator::empty() const { return _remaining_in_slice <= 0; }


void cramIterator::set(int64_t reads) {

    _curr  = 0;
    _reads = static_cast<int32_t>(reads);

}


std::shared_ptr<Codec> cramIterator::at(ExtData_t type) {

    return _container.slices[_slice].codecs.map.at(type);

}


void cramIterator::_next_slice() {

    _slice++;

    if (_slice >= _container.slices.size()) {
        _container = cramContainer(_file);
        _slice = 0;
    }

    _remaining_in_slice = _container.slices[_slice].records;

    int32_t _tmp = _container.slices[_slice].start;
    _offset      = std::max(_tmp - 1, 0);

}


void cramIterator::_next_op() {

    // The FP file tells us the delta in the query position of the
    // feature, relative to the previous feature. It uses one-based indexing
    // (i.e. the first position in the read has jump = 1) which we need to
    // account for.

    int32_t jump = at(ExtData_t::FP)->integer();
    if (_cigar.empty()) { jump--; }

    // Since matches are not stored in CRAM files, we need to determine if there
    // have been any matches since the previous read feature.

    // We can compute this by noting that the jump is given by the number
    // of matches plus the length of the query sequence consumed by the
    // previous op.

    int32_t matches = jump - _prev_qlength;
    if (matches > 0) {
        CIGAR_op match(CIGAR_t::MATCH, matches);
        _cigar.extend(match);
        _qpos += matches;
    }

    // The FC field tells us the read feature type, which we use to generate
    // the corresponding CIGAR op

    uint8_t fc     = at(ExtData_t::FC)->byte();
    ExtData_t type = _fc_to_enum(fc);
    CIGAR_op op    = _op_from_fc_enum(type, *this);
    _cigar.extend(op);

    _prev_qlength = op.qlength();
    _qpos += _prev_qlength;
    _remaining--;

}


void cramIterator::to(int32_t reference) {

    if (_reads == 0 || !_container.multi()) { return; }

    _tmp = next();
    while (_tmp.reference != reference && !end()) { _tmp = next(); _curr--; }
    _use_tmp = true;

}


std::vector<base_t> cramIterator::sub(base_t code) {

    return _container.substitution().parse(code);

}


Alignment cramIterator::next() {

    if (_use_tmp) { _use_tmp = false; return _tmp; }

    if (_remaining_in_slice <= 0) { _next_slice(); }
    _remaining_in_slice--;

    int32_t reference = UNALIGNED;
    if (_container.multi()) {
        reference = at(ExtData_t::RI)->integer();
    } else {
        reference = _container.reference;
    }

    if (reference == UNALIGNED) { return {}; }
    _curr++;

    bool aligned = true;

    // The FN field tells us the number of read features. In the case that there
    // is only one query, the FN field is not present, so we use the length of the
    // FC stream instead.

    _remaining = at(ExtData_t::FN)->integer();

    // The AP field tells us the (delta of the) alignment offset. In the case that
    // there is only one query, the AP field is not present, and we defer to the 
    // slice header's value.

    if (_container.delta()) {
        _offset += at(ExtData_t::AP)->integer();
    } else {
        _offset  = at(ExtData_t::AP)->integer() - 1;
    }

    _cigar        = CIGAR(CIGAR_RESERVE);
    _qpos         = 0;
    _prev_qlength = 0;

    while (_remaining > 0) { _next_op(); }

    // If there is a match at the end of the read, it is not accounted for
    // by jumps in the read features, so we need to access the read length
    // to determine whether it exists.

    // The read length is stored in the RL field. If there is only one
    // query, the RL field is not present, and we use the length of
    // the QS field instead.

    int32_t length = at(ExtData_t::RL)->integer();
    int32_t matches = length - _qpos;

    if (matches > 0) {
        CIGAR_op match(CIGAR_t::MATCH, matches);
        _cigar.extend(match);
    }

    // The MQ field tells us the mapping quality of the read

    qual_t mapq = at(ExtData_t::MQ)->integer();

    // The QS field tells us the PHRED scores of the read

    std::vector<uint8_t> scores = at(ExtData_t::QS)->array(length);
    PHRED phred(scores);

    return {aligned, mapq, length, _offset, reference, _cigar, phred};

}





//
// cramFile
//





static inline void _build_cram_index(BGZF* file, const std::string& name, int32_t references) {

    _throw_if_exists(name);
    std::ofstream out(name);

    IndexBlock block;
    int64_t unaligned = 0;

    block.tell(file);
    block.write_ptr(out);

    cramContainer container(file);

    int32_t tid  = 0;
    int32_t curr = 0;

    while (!container.eof()) {

        curr = container.reference;

        if (container.multi()) {

            for (auto& slice : container.slices) {

                for (int32_t ix = 0; ix < slice.records; ix++) {

                    curr = slice.index();

                    if (curr > tid) {

                        block.write_reads(out);
                        block.reads = 0;

                        for (int32_t jx = tid + 1; jx < curr; jx++) {
                            block.write_ptr(out);
                            block.write_reads(out);
                        }

                        block.reads++;
                        tid = curr;
                        block.write_ptr(out);

                    } else if (curr == UNALIGNED) {

                        unaligned++;

                    } else {

                        block.reads++;

                    }

                }

            }

        } else if (curr > tid) {

            block.write_reads(out);
            block.reads = 0;

            for (int32_t jx = tid + 1; jx < curr; jx++) {
                block.write_ptr(out);
                block.write_reads(out);
            }

            block.reads += container.aligned();
            unaligned += container.unaligned();

            tid = curr;
            block.write_ptr(out);


        } else {

            block.reads += container.aligned();
            unaligned += container.unaligned();

        }


        block.tell(file);
        container = cramContainer(file);

    }

    block.write_reads(out);

    // Write the unaligned read count

    out.write(reinterpret_cast<char*>(&unaligned), sizeof(int64_t));
    out.close();

}


cramFile::cramFile(const std::string& name) 
    : File(name, FileType::CRAM) {

    cramHeader _data = _read_cram_header(_hts_bgzf);
    _sorted     = _data.sorted;
    _references = _data.references;

    std::string _index_name = _name + CMUTS_INDEX;
    if (!std::filesystem::exists(_index_name)) {
        _build_cram_index(_hts_bgzf, _index_name, _references);
    }
    _index = Index(_index_name, _references);

}


std::shared_ptr<Iterator> cramFile::get(int32_t ix, bool seek) {

    IndexBlock block = _index.read(ix);

    // Seek to the correct alignment if necessary

    if (seek || _iterator->empty()) {

        _seek_bgzf(_hts_bgzf, block.ptr);

        _iterator = std::make_shared<cramIterator>(_hts_bgzf, block.reads);
        _iterator->to(ix);

    } 

    // Otherwise, just tell the iterator the number of reads in the alignment

    else {

        _iterator->set(block.reads);

    }

    return _iterator;

}





//
// From common.hpp
//





std::unique_ptr<File> _get_cram(const std::string& name) {

    return std::make_unique<cramFile>(name);

}





} // namespace HTS
