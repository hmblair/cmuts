#include "hdf5.hpp"

namespace HDF5 {

namespace H5 = HighFive;

static inline void __disable_logging() { H5Eset_auto2(H5E_DEFAULT, NULL, NULL); }

Manager::Manager() { __disable_logging(); }

Manager::~Manager() { H5close(); }

static inline void _throw_if_bad_compression(int compression) {
    if (compression < 0 || compression > 9) {
        throw std::runtime_error("Invalid compression level of " + std::to_string(compression) + ": must be between 0 and 9 (inclusive)");
    }
}

static inline std::string _path(const std::string& name) {
    std::filesystem::path path(name);
    return path.parent_path().string() + "/" + path.stem().string();
}

void _throw_if_object_exists(const File& file, const std::string& name) {
    std::string path = _path(name);
    if (file.exist(path)) {
        throw std::runtime_error("The object \"" + path + "\" already exists in the file \"" + file.name() + "\".");
    }
}

void _throw_if_object_exists(const File& file, const std::vector<std::string>& names) {
    for (const auto& name : names) {
        _throw_if_object_exists(file, name);
    }
}

static inline H5::FileAccessProps _create_fapl(const MPI::Manager& mpi) {

    H5::FileAccessProps fapl;

    #ifdef MPI_BUILD
    fapl.add(H5::MPIOFileAccess{mpi.comm(), mpi.info()});
    fapl.add(H5::MPIOCollectiveMetadata{});
    #endif

    return fapl;

}

File::File(
    const std::string& filename,
    AccessMode flags,
    const MPI::Manager& mpi,
    int chunk_size,
    int compression
) : H5::File(filename, flags, _create_fapl(mpi)),
    _name(filename),
    _chunk_size(chunk_size),
    _compression(compression) {

    _throw_if_bad_compression(compression);

}

std::string File::name() const {
    return _name;
}

int File::chunk_size() const {
    return _chunk_size;
}

int File::compression() const {
    return _compression;
}

template <typename dtype, size_t N>
Writer<dtype, N> File::writer(
    const std::vector<size_t>& dims,
    const std::string& name
) {
    return Writer<dtype, N>(*this, name, dims, _chunk_size, _compression);
}

template <typename dtype, size_t N>
Memspace<dtype, N> File::memspace(
    const std::vector<size_t>& dims,
    const std::string& name
) {
    return Memspace<dtype, N>(*this, dims, name);
}

static inline H5::DataSetCreateProps __get_dsprops(
    std::vector<size_t> chunkdims,
    int compression
) {

    std::vector<hsize_t> _chunkdims(chunkdims.begin(), chunkdims.end());

    H5::DataSetCreateProps dsprops;
    dsprops.add(H5::Chunking(_chunkdims));
    dsprops.add(H5::Deflate(compression));

    return dsprops;

}

static inline H5::PropertyList<(H5::PropertyType)4> __get_wprops() {

    H5::PropertyList<(H5::PropertyType)4> wprops = H5::DataTransferProps{};

    #ifdef MPI_BUILD
    wprops.add(H5::UseCollectiveIO{});
    #endif

    return wprops;

}





//
// Writer
//





template <typename dtype, size_t N>
Writer<dtype, N>::Writer(
    File& file,
    const std::string& name,
    const std::vector<size_t>& dims,
    size_t _chunksize,
    int compression
) : dims(dims), offsetdims(dims.size(), 0), _chunksize(_chunksize), wprops(__get_wprops()) {

    if (dims.size() != N) {
        throw std::runtime_error("The number of dimensions does not match the template parameter.");
    }

    chunkdims = dims;
    chunkdims[0] = std::min(static_cast<size_t>(_chunksize), dims[0]);

    auto dsprops = __get_dsprops(chunkdims, compression);
    dset = file.createDataSet<dtype>(name, H5::DataSpace(dims), dsprops);

}


template <typename dtype, size_t N>
int64_t Writer<dtype, N>::size() const { return _chunksize; }


template <typename dtype, size_t N>
void Writer<dtype, N>::write(
    const arr_t<dtype, N>& data,
    int64_t offset
) const {

    std::vector<size_t> _offsetdims = offsetdims;
    _offsetdims[0] = offset;

    dset.select(offsetdims, chunkdims)
        .write(data, wprops);

}


template <typename dtype, size_t N>
void Writer<dtype, N>::write(
    const arr_t<dtype, N>& data,
    int64_t offset,
    int64_t size
) const {

    std::vector<size_t> _offsetdims = offsetdims;
    _offsetdims[0] = offset;

    std::vector<size_t> _chunkdims = chunkdims;
    _chunkdims[0] = size;

    dset.select(_offsetdims, _chunkdims)
        .write(data, wprops);

}


template <typename dtype, size_t N>
void Writer<dtype, N>::safe_write(
    const arr_t<dtype, N>& data,
    int64_t offset
) const {

    int64_t _dim = static_cast<int64_t>(dims[0]);
    int64_t _offset = std::min(offset, _dim);
    int64_t _size = std::min(_chunksize, _dim - _offset);

    write(data, _offset, _size);

}





//
// Memspace
//





template <typename dtype, size_t N>
Memspace<dtype, N>::Memspace(
    File& file,
    const std::vector<size_t>& dims,
    const std::string& name
) : writer(file.writer<dtype, N>(dims, name)) {

    if (dims.size() != N) {
        throw std::runtime_error("The number of dimensions does not match the template parameter.");
    }

    std::vector<size_t> _buffer_dims = dims;
    _buffer_dims[0] = file.chunk_size();

    _data.resize(_buffer_dims);
    clear();

}


template <typename dtype, size_t N>
view_t<dtype, N> Memspace<dtype, N>::view(int64_t ix) {
    return xt::view(_data, ix);
}


template <typename dtype, size_t N>
void Memspace<dtype, N>::write(int64_t ix) const {
    writer.write(_data, ix);
}


template <typename dtype, size_t N>
void Memspace<dtype, N>::safe_write(int64_t ix) const {
    writer.safe_write(_data, ix);
}


template <typename dtype, size_t N>
void Memspace<dtype, N>::clear() {
    _data.fill(0.0);
}


template <typename dtype, size_t N>
int64_t Memspace<dtype, N>::size() const {
    return writer.size();
}


template <typename dtype, size_t N>
void Memspace<dtype, N>::resize(const std::vector<size_t> dims) {

    _data.resize(dims);
    clear();

}





//
// Template instantiation
//





template class Writer<float, 5>;
template class Memspace<float, 5>;
template Writer<float, 5> File::writer<float, 5>(
    const std::vector<size_t>& dims,
    const std::string& name
);
template Memspace<float, 5> File::memspace<float, 5>(
    const std::vector<size_t>& dims,
    const std::string& name
);

template class Writer<float, 4>;
template class Memspace<float, 4>;
template Writer<float, 4> File::writer<float, 4>(
    const std::vector<size_t>& dims,
    const std::string& name
);
template Memspace<float, 4> File::memspace<float, 4>(
    const std::vector<size_t>& dims,
    const std::string& name
);

template class Writer<float, 2>;
template class Memspace<float, 2>;
template Writer<float, 2> File::writer<float, 2>(
    const std::vector<size_t>& dims,
    const std::string& name
);
template Memspace<float, 2> File::memspace<float, 2>(
    const std::vector<size_t>& dims,
    const std::string& name
);

template class Writer<float, 3>;
template class Memspace<float, 3>;
template Writer<float, 3> File::writer<float, 3>(
    const std::vector<size_t>& dims,
    const std::string& name
);
template Memspace<float, 3> File::memspace<float, 3>(
    const std::vector<size_t>& dims,
    const std::string& name
);

template class Writer<int8_t, 2>;
template class Memspace<int8_t, 2>;
template Writer<int8_t, 2> File::writer<int8_t, 2>(
    const std::vector<size_t>& dims,
    const std::string& name
);
template Memspace<int8_t, 2> File::memspace<int8_t, 2>(
    const std::vector<size_t>& dims,
    const std::string& name
);


} // namespace HDF5
