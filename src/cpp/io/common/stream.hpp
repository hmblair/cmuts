#ifndef _CMUTS_IO_COMMON_STREAM_HEADER
#define _CMUTS_IO_COMMON_STREAM_HEADER

#include <cstdint>
#include <vector>
#include <string>
#include <span>

extern "C" {
    #include <htslib/bgzf.h>
    #include <zlib.h>
}

#include "types.hpp"


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

    virtual void skip(int32_t length) {};
    virtual uint8_t byte()                               = 0;
    virtual std::vector<uint8_t> bytes(int32_t length)   = 0;
    virtual std::vector<uint8_t> line(uint8_t delimiter) = 0;
    virtual std::string str(uint8_t delimiter = '\n')    = 0;

    int32_t size() const;
    int32_t remaining() const;
    bool end() const noexcept;

    int32_t bits(int32_t length);
    void memcpy(uint8_t* dest, int32_t n);
    int32_t itf8();
    int64_t ltf8();

};


//
// Streams from in-memory data
//


class DataStream : public ByteStream {
private:

    std::span<const uint8_t> _data;
    int32_t _pos = 0;

public:

    explicit DataStream(std::span<const uint8_t> data);

    uint8_t byte() override;
    std::vector<uint8_t> bytes(int32_t length) override;
    std::vector<uint8_t> line(uint8_t delimiter) override;
    std::string str(uint8_t delimiter) override;


};


class RansStream : public DataStream {
private:

    std::vector<uint8_t> _rans_data;
    explicit RansStream(std::vector<uint8_t> data);

public:

    explicit RansStream(std::span<const uint8_t> data, int32_t raw);

};


class ZlibStream : public ByteStream {
protected:

    z_stream _stream;
    std::vector<uint8_t> _buffer;

    int32_t _buffer_pos = 0;
    int32_t _buffer_end = 0;
    bool _open          = false;

public:

    ZlibStream(std::span<const uint8_t> data, int32_t raw, int32_t buffer);
    ~ZlibStream() override;

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


class BgzfFileStream : public ByteStream {
private:

    BGZF* _file;

public:

    BgzfFileStream(BGZF* file, int32_t size);

    void skip(int32_t length) override;
    uint8_t byte() override;
    std::vector<uint8_t> bytes(int32_t length) override;
    std::vector<uint8_t> line(uint8_t delimiter) override;
    std::string str(uint8_t delimiter) override;

};


class ZlibFileStream : public ZlibStream {
private:

    BgzfFileStream _bgzf;
    std::vector<uint8_t> _bgzf_buffer;
    std::span<const uint8_t> _buffer_span;

    int32_t _bgzf_buffer_pos = 0;
    int32_t _in_remaining    = 0;

public:

    ZlibFileStream(BGZF* file, std::span<const uint8_t> data, int32_t size, int32_t raw, int32_t buffer);
    ZlibFileStream(BGZF* file, int32_t size, int32_t raw, int32_t buffer);

    void fill() override;

};


#endif
