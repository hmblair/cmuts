#include "h5utils.h"

Filespace makeFilespace(const char* filename) {

    Filespace filespace;
    // By specifying H5F_ACC_EXCL, an error will be thrown if the file
    // already exists.
    filespace.fileID = H5Fcreate(filename, H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
    if (filespace.fileID < 0) {
        fprintf(stderr, "Error creating HDF5 file\n");
        exit(EXIT_FAILURE);
    }

    filespace.filename = filename;

    return filespace;

}


Filespace openFilespace(const char* filename, unsigned flags) {

    Filespace filespace;

    filespace.fileID = H5Fopen(filename, flags, H5P_DEFAULT);
    if (filespace.fileID < 0) {
        fprintf(stderr, "Error opening file %s\n", filename);
        exit(EXIT_FAILURE);
    }

    filespace.filename = filename;

    return filespace;

}


void closeFilespace(Filespace filespace) {

    herr_t status = H5Fclose(filespace.fileID);
    if (status < 0) {
        fprintf(stderr, "Error closing the file %s.\n", filespace.filename);
        exit(EXIT_FAILURE);
    }

}


bool groupExists(Filespace filespace, const char *name) {

    htri_t status = H5Lexists(filespace.fileID, name, H5P_DEFAULT);
    // The group exists
    if (status > 0) {
        return true;
    } 
    // The group does not exist
    else if (status == 0) {
        return false;
    }
    // There was an error
    else {
        fprintf(stderr, "There was an error checking if the group \"%s\" exists.\n", name);
        exit(EXIT_FAILURE);
    }

}


Dataspace makeDataspace(
    Filespace filespace, 
    const char *name, 
    size_t ndims, 
    hsize_t *dims, 
    DataspaceFlags dFlags
) {

    if (dFlags.compression < 0 || dFlags.compression > 9) {
        fprintf(stderr, "The compression level (%d) must be between 0 and 9 inclusive.\n", dFlags.compression);
        exit(EXIT_FAILURE);
    }

    Dataspace dataspace;
    dataspace.ndims = ndims;
    // Copy the dimensions into their own memory space.
    dataspace.dims = calloc(ndims, sizeof(hsize_t));
    if (!dataspace.dims) {
        fprintf(stderr, "There was an error allocating memory.\n");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < ndims; i++) {
        dataspace.dims[i] = dims[i];
    }
    dataspace.name = name;

    // By default, the location of the dataset will be in root.
    hid_t locationID = filespace.fileID;
    // If a group name is specified, then check if the group 
    // exists in the file. If it does, open it, and check if
    // there is already a dataset with this name. Else, create it.
    if (dFlags.group) {
        if (groupExists(filespace, dFlags.group)) {
            locationID = H5Gopen(filespace.fileID, dFlags.group, H5P_DEFAULT);
            if (locationID < 0) {
                fprintf(stderr, "There was an error opening the group \"%s\".\n", dFlags.group);
                exit(EXIT_FAILURE);
            }
        } else {
            locationID = H5Gcreate(filespace.fileID, dFlags.group, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            if (locationID < 0) {
                fprintf(stderr, "There was an error creating the group \"%s\".\n", dFlags.group);
                exit(EXIT_FAILURE);
            }
        }
    }

    hid_t plist = H5Pcreate(H5P_DATASET_CREATE);
    if (plist < 0) {
        fprintf(stderr, "Failed to create dataset property list for %s\n", name);
        exit(EXIT_FAILURE);
    }
    // Set chunking, in order to use compression
    hsize_t *chunk = calloc(ndims, sizeof(hsize_t));
    if (!chunk) {
        fprintf(stderr, "There was an error allocating memory.\n");
        exit(EXIT_FAILURE);
    }
    hsize_t *maxDims = calloc(ndims, sizeof(hsize_t));
    if (!maxDims) {
        fprintf(stderr, "There was an error allocating memory.\n");
        exit(EXIT_FAILURE);
    }
    chunk[0] = dFlags.chunkSize;
    maxDims[0] = H5S_UNLIMITED;
    for (int i = 1; i < ndims; i++) {
        chunk[i] = dims[i];
        maxDims[i] = dims[i];
    }
    herr_t chunkStatus = H5Pset_chunk(plist, ndims, chunk);
    if (chunkStatus < 0) {
        fprintf(stderr, "Failed to set the chunk size of %d in %s.\n", dFlags.chunkSize, name);
        exit(EXIT_FAILURE);
    }
    H5Pset_deflate(plist, dFlags.compression);
    H5Pset_shuffle(plist);
    H5Pset_fletcher32(plist);

    #pragma omp critical 
    {

        // Create a dataspace to store the data
        dataspace.dataspaceID = H5Screate_simple(ndims, dims, maxDims);
        if (dataspace.dataspaceID< 0) {
            fprintf(stderr, "Error creating dataspace %s in file %s\n", name, filespace.filename);
            exit(EXIT_FAILURE);
        }
        // Create the corresponding dataset
        dataspace.datasetID = H5Dcreate2(locationID, name, H5T_NATIVE_FLOAT, dataspace.dataspaceID, H5P_DEFAULT, plist, H5P_DEFAULT);
        if (dataspace.datasetID< 0) {
            fprintf(stderr, "Error creating dataset %s in file %s\n", name, filespace.filename);
            exit(EXIT_FAILURE);
        }

    }

    H5Pclose(plist);
    free(chunk);
    free(maxDims);
    // If a group was opened, close it after creating the dataset
    if (locationID != filespace.fileID) {
        H5Gclose(locationID);
    }

    return dataspace;

}


Dataspace openDataspace(Filespace filespace, const char* name) {

    Dataspace dataspace;
    dataspace.name = name;
    dataspace.datasetID = H5Dopen2(filespace.fileID, name, H5P_DEFAULT);
    if (dataspace.datasetID < 0) {
        fprintf(stderr, "Error opening dataset %s\n", name);
        exit(EXIT_FAILURE);
    }

    // Get the dataspace of the dataset
    dataspace.dataspaceID = H5Dget_space(dataspace.datasetID);

    // Get the number of dimensions in the dataset
    dataspace.ndims = H5Sget_simple_extent_ndims(dataspace.dataspaceID);
    if (dataspace.ndims < 1) {
        fprintf(stderr, "Dataset has no dimensions.\n");
        exit(EXIT_FAILURE);
    }

    dataspace.dims = calloc(dataspace.ndims, sizeof(hsize_t));
    if (!dataspace.dims) {
        fprintf(stderr, "There was an error allocating memory.\n");
        exit(EXIT_FAILURE);
    }
    H5Sget_simple_extent_dims(dataspace.dataspaceID, dataspace.dims, NULL);

    return dataspace;

}


void closeDataspace(Dataspace dataspace) {

    free(dataspace.dims);
    herr_t status = H5Sclose(dataspace.dataspaceID);
    if (status < 0) {
        fprintf(stderr, "Error closing the dataspace in %s.\n", dataspace.name);
        exit(EXIT_FAILURE);
    }
    status = H5Dclose(dataspace.datasetID);
    if (status < 0) {
        fprintf(stderr, "Error closing the dataset in %s.\n", dataspace.name);
        exit(EXIT_FAILURE);
    }

}


DataspaceFlags defaultDataspaceFlags() {
    DataspaceFlags dFlags;
    dFlags.chunkSize = CHUNK_SIZE;
    dFlags.compression = DEFAULT_COMPRESSION;
    dFlags.group = NULL;
    return dFlags;
}


Memspace initMemspace(Dataspace dataspace, size_t slabSize) {

    if (slabSize > dataspace.dims[0]) {
        fprintf(stderr, "The chosen slice size (%lu) is larger than the size of the dataset (%llu).\n", slabSize, dataspace.dims[0]);
        exit(EXIT_FAILURE);
    }

    Memspace memspace;
    memspace.ndims = dataspace.ndims;
    memspace.dataSize = sizeof(float);
    // Initialise the slice to be pointing at the beginning of
    // the dataset
    memspace.slabOffsets = calloc(dataspace.ndims, sizeof(hsize_t));
    if (!memspace.slabOffsets) {
        fprintf(stderr, "There was an error allocating memory.\n");
        exit(EXIT_FAILURE);
    }
    // The slice has a size of slabSize along the first
    // dimension and a full extend along all the others
    memspace.slabDims = calloc(dataspace.ndims, sizeof(hsize_t));
    if (!memspace.slabDims) {
        fprintf(stderr, "There was an error allocating memory.\n");
        exit(EXIT_FAILURE);
    }
    memspace.slabDims[0] = slabSize;
    for (int i = 1; i < dataspace.ndims; i++) {
        memspace.slabDims[i] = dataspace.dims[i];
    }
    // Allocate memory for one slice. To do so, we need to get the
    // total number of elements
    memspace.numElements = slabSize;
    for (int i = 1; i < dataspace.ndims; ++i) {
        memspace.numElements *= dataspace.dims[i];
    }
    memspace.data = (float*)calloc(memspace.numElements, sizeof(float));
    if (!memspace.data) {
        fprintf(stderr, "Memory allocation failed when initialising the memspace.\n");
        exit(EXIT_FAILURE);
    }

    #pragma omp critical 
    {

        // Define a memory space the same size as the slice
        memspace.memspaceID = H5Screate_simple(memspace.ndims, memspace.slabDims, NULL);
        if (memspace.memspaceID < 0) {
            fprintf(stderr, "Error creating memspace for the dataspace %s\n", dataspace.name);
            exit(EXIT_FAILURE);
        }

    }

    return memspace;

}

void resizeMemspace(Memspace *memspace, size_t slabSize) {


    if (slabSize == memspace->slabDims[0]) {
        return;
    }

    memspace->slabDims[0] = slabSize;

    memspace->numElements = 1;
    for (int i = 0; i < memspace->ndims; ++i) {
        memspace->numElements *= memspace->slabDims[i];
    }
    memspace->data = realloc(memspace->data, memspace->numElements * sizeof(float));
    if (memspace->data == NULL) {
        fprintf(stderr, "Error reallocating memory when resizing a memspace.\n");
        exit(EXIT_FAILURE);
    }
    hsize_t *memspaceStart = calloc(memspace->ndims, sizeof(hsize_t));

    #pragma omp critical 
    {
        // Select the hyperslab in the memory's dataspace (for consistency)
        herr_t status = H5Sselect_hyperslab(memspace->memspaceID, H5S_SELECT_SET, memspaceStart, NULL, memspace->slabDims, NULL);
        if (status < 0) {
            fprintf(stderr, "Error selecting hyperslab in memory space\n");
            exit(EXIT_FAILURE);
        }
    }

    free(memspaceStart);


}


void fillMemspace(Memspace memspace, float value) {
    
    for (int i = 0; i < memspace.numElements; i++) {
        memspace.data[i] = value;
    }

}


void freeMemspace(Memspace memspace) {

    free(memspace.data);
    free(memspace.slabDims);
    free(memspace.slabOffsets);

    #pragma omp critical 
    {
        int status = H5Sclose(memspace.memspaceID);
        if (status < 0) {
            fprintf(stderr, "There was an error closing the memspace.\n");
            exit(EXIT_FAILURE);
        }
    }

}

// Load the slab starting from index ix into the memory
// allocated in the memspace
void loadSlab(Dataspace dataspace, Memspace memspace, hsize_t ix) {

    hsize_t maxIX = dataspace.dims[0] - memspace.slabDims[0]; 
    if (ix > maxIX) {
        fprintf(stderr, "Index out of bounds: %llu (max: %llu)\n", ix, maxIX);
        exit(EXIT_FAILURE);
    }

    // Move the memspace offset to that it points to the beginning
    // of the desired slice
    memspace.slabOffsets[0] = ix;

    // Select the hyperslab
    H5Sselect_hyperslab(dataspace.dataspaceID, H5S_SELECT_SET, memspace.slabOffsets, NULL, memspace.slabDims, NULL);

    // Read the selected record into memory
    herr_t status = H5Dread(dataspace.datasetID, H5T_NATIVE_FLOAT, memspace.memspaceID, dataspace.dataspaceID, H5P_DEFAULT, memspace.data);
    if (status < 0) {
        fprintf(stderr, "Error reading slab for index %llu\n", ix);
        exit(EXIT_FAILURE);
    }

}


void writeSlab(Dataspace dataspace, Memspace memspace, hsize_t ix) {

    hsize_t maxIX = dataspace.dims[0] - memspace.slabDims[0]; 
    if (ix > maxIX) {
        fprintf(stderr, "Index out of bounds: %llu (max: %llu)\n", ix, maxIX);
        exit(EXIT_FAILURE);
    }

    // Move the memspace offset so that it points to the beginning
    // of the desired slice
    memspace.slabOffsets[0] = ix;


    #pragma omp critical 
    {

        // Select the hyperslab in the file's dataspace
        herr_t status = H5Sselect_hyperslab(dataspace.dataspaceID, H5S_SELECT_SET, memspace.slabOffsets, NULL, memspace.slabDims, NULL);
        if (status < 0) {
            fprintf(stderr, "Error selecting hyperslab for index %llu\n", ix);
            exit(EXIT_FAILURE);
        }

        // Write the data from the memspace to the selected hyperslab in the file
        status = H5Dwrite(dataspace.datasetID, H5T_NATIVE_FLOAT, memspace.memspaceID, dataspace.dataspaceID, H5P_DEFAULT, memspace.data);
        if (status < 0) {
            fprintf(stderr, "Error writing slab for index %llu\n", ix);
            exit(EXIT_FAILURE);
        }

    }

}


bool isFloat(Dataspace dataspace) {

    hid_t dtype = H5Dget_type(dataspace.datasetID);

    if (H5Tget_class(dtype) != H5T_FLOAT) {
        printf("%d\n",H5Tget_class(dtype));
        H5Tclose(dtype);
        return false;
    }

    return true;
}
