#ifndef HDF5_HEADER
#define HDF5_HEADER

#include <ranges>
#include <stdexcept>
#include <string>
#include <vector>
#include <span>

#include <highfive/highfive.hpp>
#include <highfive/xtensor.hpp>
#include <xtensor/containers/xarray.hpp>
#include <xtensor/io/xio.hpp>
#include <xtensor/views/xview.hpp>

#include "infra/mpi.hpp"
#include "infra/utils.hpp"

// The array type used internally
template <typename dtype, size_t N>
using arr_t = xt::xtensor<dtype, N>;
// Most arrays are passed as views
template <typename dtype, size_t N>
using view_t = decltype(xt::view(std::declval<arr_t<dtype, N>&>(), std::declval<int32_t>()));

namespace HDF5 {

namespace H5 = HighFive;

const H5::File::AccessMode R   = H5::File::ReadOnly;
const H5::File::AccessMode RW  = H5::File::ReadWrite;
const H5::File::AccessMode RWC = H5::File::ReadWrite | H5::File::Create;

class Manager{
public:

    Manager();
    ~Manager();

};

template <typename dtype, size_t N>
class Writer;

template <typename dtype, size_t N>
class Memspace;

class File : public H5::File {
private:

    std::string _name;
    int _chunk_size  = 0;
    int _compression = 0;

public:

    File(
        const std::string& filename,
        AccessMode flags,
        const MPI::Manager& mpi,
        int chunk_size,
        int compression
    );

    std::string name() const;
    int chunk_size() const;
    int compression() const;

    template <typename dtype, size_t N>
    Writer<dtype, N> writer(
        const std::vector<size_t>& dims,
        const std::string& name
    );

    template <typename dtype, size_t N>
    Memspace<dtype, N> memspace(
        const std::vector<size_t>& dims,
        const std::string& name
    );

};

template <typename dtype, size_t N>
class Writer {
private:

    H5::DataSet dset;
    std::vector<size_t> dims;
    std::vector<size_t> chunkdims;
    std::vector<size_t> offsetdims;
    int32_t _chunksize;
    H5::PropertyList<(H5::PropertyType)4> wprops;

public:

    Writer(
        File& file,
        const std::string& name,
        const std::vector<size_t>& dims,
        size_t _chunksize,
        int compression
    );

    int32_t size() const;

    void write(
        const arr_t<dtype, N>& data,
        int32_t offset
    ) const;

    void write(
        const arr_t<dtype, N>& data,
        int32_t offset,
        int32_t size
    ) const;

    void safe_write(
        const arr_t<dtype, N>& data,
        int32_t offset
    ) const;


};

template <typename dtype, size_t N>
class Memspace {
private:

    Writer<dtype, N> writer;
    arr_t<dtype, N> _data;

public:

    Memspace(
        File& h5file,
        const std::vector<size_t>& dims,
        const std::string& name
    );

    int32_t size() const;
    view_t<dtype, N> view(int32_t ix);
    void write(int32_t ix) const;
    void safe_write(int32_t ix) const;
    void clear();
    void resize(const std::vector<size_t> dims);

};

void _throw_if_object_exists(const File& file, const std::string& name);
void _throw_if_object_exists(const File& file, const std::vector<std::string>& names);

} // namespace HDF5

#endif
