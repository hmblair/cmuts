#include "main.hpp"

static inline void __write_sequences(
    const HTS::FASTA& fasta,
    HDF5::File& hdf5,
    const MPI::Manager& mpi
) {

    if (!hdf5.exist(SEQUENCE_DS)) {
        try {
            if (mpi.root()) {
                std::cout << "Tokenizing " << fasta.size() << " sequences.\n";
            }
            fasta.to_hdf5(hdf5, mpi);
        } catch (const std::exception& e) {
            mpi.err() << "Error: " << e.what() << "\n";
        }
    } else {
        mpi.err() << "Warning: The file \"" << hdf5.name() << "\" already contains the dataset \"" << SEQUENCE_DS << "\". No tokenization can be done.\n";
    }

}

static inline void __cleanup(
    const MPI::Manager& mpi,
    const cmutsProgram& opt
) {
    try {
        if (opt.overwrite) {
            mpi.remove(opt.output);
        }
    } catch (const std::exception& e) {
        mpi.err() << "Error during cleanup: " << e.what() << "\n";
    }
}

int main(int argc, char** argv) {

    // Initialize MPI processes
    MPI::Manager mpi(argc, argv);
    // Initialize HDF5 manager
    HDF5::Manager _hdf5_manager;
    // Disable native HTS logging as it does not work well in the
    // multi-threaded environment
    HTS::__disable_logging();
    // For printing integers with commas
    __imbue();

    // Parse command line arguments
    cmutsProgram opt;
    try {
        opt.parse(argc, argv);
    } catch (const std::exception& err) {
        mpi.err() << "Error: " << err.what() << "\n";
        return EXIT_FAILURE;
    }

    // Delete the exiting output file if specified
    try {
        if (opt.overwrite) {
            mpi.remove(opt.output);
        }
    } catch (const std::exception& e) {
        mpi.err() << "Error: " << e.what() << "\n";
        return EXIT_FAILURE;
    }

    // Create the HDF5 file
    std::optional<HDF5::File> _hdf5;
    try {
        _hdf5.emplace(opt.output, HDF5::RWC, mpi, opt.chunk_size, opt.compression);
    } catch (const std::exception& e) {
        mpi.err() << "Error: " << e.what() << "\n";
        return EXIT_FAILURE;
    }
    HDF5::File hdf5 = std::move(_hdf5.value());

    try {
        HDF5::_throw_if_object_exists(hdf5, opt.files);
    } catch (const std::exception& e) {
        mpi.err() << "Error: " << e.what() << "\n";
        return EXIT_FAILURE;
    }

    // Open the reference FASTA
    std::optional<HTS::FASTA> _fasta;
    try {
        _fasta.emplace(opt.fasta);
    } catch (const std::exception& e) {
        mpi.err() << "Error: " << e.what() << "\n";
        __cleanup(mpi, opt);
        return EXIT_FAILURE;
    }
    HTS::FASTA fasta = std::move(_fasta.value());

    size_t total = opt.files.value().size();
    if (total == 0) {
        __write_sequences(fasta, hdf5, mpi);
        return EXIT_SUCCESS;
    }

    // Open the HTS files
    std::optional<HTS::FileGroup> _files;
    try {
        _files.emplace(opt.files, fasta, mpi);
    } catch (const std::exception& e) {
        mpi.err() << "Error: " << e.what() << "\n";
        __cleanup(mpi, opt);
        return EXIT_FAILURE;
    }
    HTS::FileGroup files = std::move(_files.value());
    mpi.barrier();

    // Get the desired operation mode
    cmuts::DetailLevel detail;
    try {
        detail = cmuts::detail(opt.fast, opt.joint);
    } catch (const std::exception& e) {
        mpi.err() << "Error: " << e.what() << "\n";
        __cleanup(mpi, opt);
        return EXIT_FAILURE;
    }

    // Initialise the parameters for the main counting function
    cmuts::Params params = {
        opt.min_mapq,
        opt.min_quality,
        opt.min_length,
        opt.max_length,
        opt.max_indel_length,
        opt.spread,
        opt.quality_window,
        detail,
    };

    // Initialise the stats tracker, and print the header
    cmuts::Stats stats(
        files.reads(),
        files.unaligned_reads(),
        files.references(),
        mpi
    );
    stats.header();

    size_t processed = 0;
    for (auto& input : files) {
        std::unique_ptr<cmuts::__Main> main;
        try {
            main = cmuts::get_main(input, hdf5, mpi, params, stats);
            main->run();
            processed++;
        } catch (const std::exception& e) {
            mpi.err() << "Error processing the file \"" << input.name() << "\": " << e.what() << "\n";
            continue;
        }
    }

    if (mpi.root()) {
        Utils::cursor_down(1);
    }

    if (processed > 0) {
        __write_sequences(fasta, hdf5, mpi);
    } else {
        __cleanup(mpi, opt);
    }

    if (processed < files.size()) {
        mpi.err() << "WARNING: only " << processed << " of " << files.size() << " files were processed.\n";
    }

    return EXIT_SUCCESS;

}
