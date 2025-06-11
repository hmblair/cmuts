#ifndef _CMUTS_FASTA_HEADER
#define _CMUTS_FASTA_HEADER

#include "common.hpp"
#include "hdf5.hpp"
#include "mpi.hpp"

// Indices of the bases in binary sequences

const uint8_t BIN_A   = 0x0;
const uint8_t BIN_C   = 0x1;
const uint8_t BIN_G   = 0x2;
const uint8_t BIN_T   = 0x3;
const uint8_t BIN_UNK = 0x4;

// Other binary sequence data constants

const int32_t BASES_PER_BYTE = 4;
const int32_t EOH = -1;

// HDF5 constants

const std::string FASTA_DATASET = "sequence";
const size_t FASTA_DATASET_SIZE = 2;





//
// FASTA
//





struct Offset {

    int32_t length = 0;
    int32_t offset = 0;

};


class HeaderBlock{
public:

    int32_t length    = 0;
    int32_t sequences = 0;
    bool empty() const;

};





class Header{
private:

    std::vector<HeaderBlock> _blocks;
    int32_t _offset = 0;

public:

    Header() = default;
    explicit Header(std::ifstream& file);

    int32_t size() const;
    int32_t length(int32_t ix) const;
    int32_t longest() const;
    Offset offset(int32_t ix);

};





class BinaryFASTA {
private:

    std::string _fasta_name;
    std::string _name;
    std::ifstream _file;
    Header _header;

public:

    explicit BinaryFASTA(const std::string& fasta);
    BinaryFASTA(BinaryFASTA&& other) noexcept;
    BinaryFASTA& operator=(BinaryFASTA&& other) noexcept;
    BinaryFASTA(const BinaryFASTA&) = delete;
    BinaryFASTA& operator=(const BinaryFASTA&) = delete;
    ~BinaryFASTA() = default;

    std::string name() const;
    int32_t size() const;

    int32_t length(int32_t ix) const;
    seq_t sequence(int32_t ix);
    int32_t longest() const;

    void hdf5(HDF5::File& hdf5, const MPI::Manager& mpi);

};





#endif
