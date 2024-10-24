#ifndef H5_HEADER
#define H5_HEADER

#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>

// HDF5
#include <hdf5.h>
// OMP
#include <mpi.h>

// HDF5 parameters
#define HDF5_SUFFIX ".h5"
#define CHUNK_SIZE 100
#define DEFAULT_COMPRESSION 3

typedef struct {
    const char *filename;
    hid_t fileID;
} Filespace;

typedef struct {
    const char* name;
    hid_t datasetID;
    hid_t dataspaceID;
    int ndims;
    hsize_t *dims;
} Dataspace;

typedef struct {
    hid_t memspaceID;
    hsize_t *slabDims;
    hsize_t *slabOffsets;
    int ndims;
    float *data;
    size_t dataSize;
    size_t numElements;
} Memspace;

typedef struct {
    int compression;
    int chunkSize;
    const char *group;
} DataspaceFlags;

Filespace makeFilespace(
    const char* filename, 
    MPI_Comm comm, 
    MPI_Info info
);

Filespace openFilespace(
    const char* filename, 
    unsigned flags, 
    MPI_Comm comm, 
    MPI_Info info
);

Filespace getFilespace(
    const char *filename,
    bool exists,
    unsigned flags,
    MPI_Comm comm,
    MPI_Info info
);

void closeFilespace(Filespace filespace);

bool groupExists(Filespace filespace, const char *name);

Dataspace makeDataspace(
    Filespace filespace,
    const char* name,
    size_t ndims,
    hsize_t *dims,
    DataspaceFlags flags
);

Dataspace openDataspace(Filespace filespace, const char* name);

void closeDataspace(Dataspace dataspace);

DataspaceFlags defaultDataspaceFlags();

Memspace initMemspace(Dataspace dataspace, size_t slabSize);

void resizeMemspace(Memspace *memspace, size_t slabSize);

void fillMemspace(Memspace memspace, float value);

void freeMemspace(Memspace memspace);

void loadSlab(
    Dataspace dataspace,
    Memspace memspace,
    hsize_t ix
);

void writeSlab(
    Dataspace dataspace,
    Memspace memspace,
    hsize_t ix
);

bool isFloat(Dataspace dataspace);

#endif
