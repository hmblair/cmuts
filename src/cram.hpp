#ifndef _CMUTS_CRAM_HEADER
#define _CMUTS_CRAM_HEADER

#include "common.hpp"
#include <unordered_map>

const int32_t FILE_ID_SIZE = 20;

const int32_t SM_SIZE  = 5;
const int32_t MD5_SIZE = 16;
const int32_t KEY_SIZE = 2;
const int32_t CRC_SIZE = 4;

const int32_t CRAM_EOF_POS  = 4542278;
const int32_t CRAM_MULTI    = -2;





//
// ITF8 and LTF8 integers
//





struct ITF8 {

    int32_t size  = 1;
    int32_t value = 0;
    int32_t read  = 1;

};


struct LTF8 {

    int32_t size  = 1;
    int64_t value = 0;
    int32_t read  = 1;

};


struct ITF8arr {

    int32_t size = 1;
    std::vector<int32_t> values;

};


struct LTF8arr {

    int32_t size = 1;
    std::vector<int64_t> values;

};






namespace HTS {





struct cramHeader : public Header {

    int32_t version;

};


enum class ExtData_t : uint8_t {

    BF,
    CF,
    RI,
    RL,
    AP,
    RG,
    RN,
    MF,
    NS,
    NP,
    TS,
    NF,
    TL,
    FN,
    FC,
    FP,
    DL,
    BB,
    QQ,
    BS,
    IN,
    RS,
    PD,
    HC,
    SC,
    MQ,
    BA,
    QS

};


// The number of external data types.

const int32_t AUX_DATA_SIZE = 28;

// Which external data to skip when reading a CRAM file.
// Currently: read names, ??

const std::unordered_set<ExtData_t> EXT_SKIP = {ExtData_t::RN, ExtData_t::TL};


enum class AuxDataEncoding : uint8_t {

    Int32,
    Byte,
    ByteArray

};


enum class Codec_t : uint8_t {

    ID,           // Just an ITF8 integer
    Huffman,      // Canonical Huffman coding
    ByteArrayLen, // An integer encoding of lengths followed by a byte encoding
    ByteStop,     // Stop character + ITF8 integer
    Beta,         // Offset + length, both ITF8
    SubExp        // Offset + exponent, both ITF8

};


enum class CompressionMethod : uint8_t {

    None,
    gzip,
    rans4x8

};


enum class cramBlockType : uint8_t {

    Header,      // File header
    Compression, // Compression header
    Slice,       // Slice header
    RESERVED,    // Unused
    External,    // External data
    Core         // Core data

};





//
// Headers
//





struct PreservationMap {

    bool RN = true;
    bool AP = true;
    bool RR = true;
    uint8_t SM[SM_SIZE];
    std::vector<uint8_t> TD;

};


class Codec {
protected:

    Codec_t _type;
    int32_t _id = -1;
    std::shared_ptr<ByteStream> _stream = nullptr;

public:

    explicit Codec(Codec_t type);
    Codec(Codec_t type, int32_t id);
    virtual ~Codec() = default;

    Codec_t type() const;
    int32_t id() const;
    void attach(const std::shared_ptr<ByteStream>& stream);

    virtual uint8_t byte();
    virtual int32_t integer();
    virtual std::vector<uint8_t> array();
    virtual std::vector<uint8_t> array(int32_t length);

};


class ExternalCodec : public Codec {
private:

public:

    explicit ExternalCodec(int32_t id);

    std::vector<uint8_t> array(int32_t length) override;

};


class ByteArrayStopCodec : public Codec {
private:

    uint8_t _stop = 0;

public:

    explicit ByteArrayStopCodec (int32_t id, uint8_t stop);

    std::vector<uint8_t> array() override;

};


class HuffmanCodec: public Codec {
private:

    std::vector<int32_t> _alphabet;
    std::vector<int32_t> _lengths;

public:

    HuffmanCodec(const std::vector<int32_t>& alphabet, const std::vector<int32_t>& lengths);

    uint8_t byte() override;
    int32_t integer() override;

};


class BetaCodec: public Codec {
private:

    int32_t _offset = 0;
    int32_t _length = 0;

public:

    BetaCodec(int32_t offset, int32_t length);

    int32_t integer() override;

};


class SubExpCodec: public Codec {
private:

    int32_t _offset = 0;
    int32_t _order  = 0;

public:

    SubExpCodec(int32_t offset, int32_t order);

    int32_t integer() override;

};


struct CodecMap {

    std::unordered_map<ExtData_t, std::shared_ptr<Codec>> map;
    std::unordered_map<int32_t, ExtData_t> externals;

};



//
// cramBlock
//


class SubstitutionMatrix {
private:

    std::vector<std::vector<base_t>> _decoder;

public:

    SubstitutionMatrix() = default;
    explicit SubstitutionMatrix(uint8_t* matrix);

    std::vector<base_t> parse(base_t code);

};


class cramBlock {
private:

    bool _has_crc;

public:

    explicit cramBlock(BGZF* file, bool crc);
    explicit cramBlock(uint8_t*& _data, bool crc);
    ~cramBlock() = default;

    // The block header as per the specification

    CompressionMethod method;
    cramBlockType type;
    int32_t content;
    int32_t size;
    int32_t raw;

    // The data contained in the block; exactly one will be non-null/empty

    std::span<const uint8_t> data;
    BGZF* file = nullptr;

    uint8_t byte();
    std::vector<uint8_t> bytes(int32_t n);

    void memcpy(uint8_t* dest, int32_t size);

    int32_t itf8();
    int64_t ltf8();
    std::vector<int32_t> itf8_arr();
    std::vector<int64_t> ltf8_arr();

    template <typename dtype>
    std::vector<dtype> array();

    int32_t length(uint8_t stop = 0);
    std::string str(int32_t length);

    // Are we at the end of the block

    bool end() const;

};


class cramSlice {
public:

    // Load just the header, or the header and all blocks

    explicit cramSlice(uint8_t*& data, bool crc);
    cramSlice(uint8_t*& data, CodecMap codecs, bool crc);

    int32_t reference;
    int32_t start;
    int32_t span;
    int32_t records;
    int64_t counter;
    int32_t nblocks;
    int32_t emb_ref_id;
    uint8_t md5[MD5_SIZE];

    // All codecs associated with the slice

    CodecMap codecs;

    int32_t index();

    // Is this a multi-reference slice

    bool multi() const;

    // The number of aligned reads in the slice

    int64_t aligned();

    // The number of unaligned reads in the slice

    int64_t unaligned();

};


class CompressionHeader {
private:

    PreservationMap _map;
    CodecMap _codecs;
    SubstitutionMatrix _subs;

public:

    CompressionHeader() = default;
    explicit CompressionHeader(uint8_t*& data, bool crc);

    CodecMap codecs() const;
    SubstitutionMatrix substitution() const;
    bool delta() const;

};





//
// cramContainerBase; cramContainer; cramHeaderContainer
//





class cramContainerBase {
protected:

    BGZF* _file   = nullptr;
    bool _has_crc = false;

public:

    cramContainerBase() = default;
    explicit cramContainerBase(BGZF* file, int32_t version);

    // The location of the container in the file
    int64_t ptr = 0;

    // The container header as per the specification

    int32_t size      = 0;
    int32_t reference = 0;
    int32_t offset    = 0;
    int32_t span      = 0;
    int32_t records   = 0;
    int64_t counter   = 0;
    int64_t bases     = 0;
    int32_t blocks    = 0;
    std::vector<int32_t> landmarks;
    int32_t crc32     = 0;

    // Is this the EOF container

    bool eof() const;

    // Is this a multi-reference container

    bool multi() const;

};


class cramContainer : public cramContainerBase {
private:

    // The compression header; only present if !eof()

    CompressionHeader _comp;

public:

    cramContainer() = default;
    explicit cramContainer(BGZF* file, int32_t version);

    // The raw data of the container

    std::vector<uint8_t> data;

    // The slices in the container

    std::vector<cramSlice> slices;
    cramSlice slice(int32_t ix);

    // The number of aligned reads in the container

    int64_t aligned();

    // The number of unaligned reads in the container

    int64_t unaligned();

    // Get all initialized codecs

    CodecMap codecs() const;

    SubstitutionMatrix substitution() const;

    bool delta() const;

};


class cramHeaderContainer : public cramContainerBase {
private:

public:

    explicit cramHeaderContainer(BGZF* file, int32_t version);

    void skip();

    // The raw data of the container

    std::unique_ptr<ByteStream> stream;

};








class cramIterator : public Iterator {
private:

    // The BGZF file and version

    BGZF* _file     = nullptr;
    int32_t _version = 0;

    // Reference to the blocks in the slice

    cramContainer _container;

    // Load the next slice, which may be in the next container

    int32_t _slice              = 0;
    int32_t _remaining_in_slice = 0;
    void _next_slice();

    // Internal parameters for constructing the CIGAR

    CIGAR _cigar;

    int32_t _qpos         = 0;
    int32_t _prev_qlength = 0;
    int32_t _offset       = 0;
    int32_t _remaining    = 0;

    // Add an op to the CIGAR

    void _next_op();

    // Used in the to() method

    Alignment _tmp;
    bool _use_tmp = false;

public:

    explicit cramIterator(BGZF* file, int32_t reads, int32_t version);
    cramIterator(cramIterator&&) = delete;
    cramIterator& operator=(cramIterator&&) = delete;

    std::shared_ptr<Codec> at(ExtData_t type);
    std::vector<base_t> sub(base_t code);
    Alignment next() override;
    void to(int32_t reference);
    bool empty() const;
    void set(int64_t reads);

};





//
// cramFile
//





class cramFile : public File {
private:

    int32_t _version = 0;
    std::shared_ptr<cramIterator> _iterator;

public:

    explicit cramFile(const std::string& name);
    ~cramFile() override = default;

    std::shared_ptr<Iterator> get(int32_t ix, bool seek) override;

};





//
// Printing overloads
//





std::ostream& operator<<(std::ostream& os, const cramBlock& block);
std::ostream& operator<<(std::ostream& os, const cramContainer& container);





} // namespace HTS




#endif
