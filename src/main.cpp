#include "main.hpp"

static inline void _print_title(const MPI::Manager& mpi) {

    mpi.down();
    if (mpi.root()) { Utils::version(); }
    mpi.divide();

}


static inline void __write_sequences(
    BinaryFASTA& fasta,
    HDF5::File& hdf5,
    const MPI::Manager& mpi
) {

    mpi.down();

    if (!hdf5.exist(FASTA_DATASET)) {
        try {
            fasta.hdf5(hdf5, mpi);
            mpi.out() << "        Tokenized " << fasta.size() << " sequences.\n";
        } catch (const std::exception& e) {
            mpi.err() << "Error: " << e.what() << "\n";
        }
    } else {
        mpi.err() << "WARNING: The file \"" << hdf5.name() << "\" already contains the dataset \"" << FASTA_DATASET << "\". No tokenization can be done.\n";
    }

}


static inline void __cleanup(
    const MPI::Manager& mpi,
    const cmutsProgram& opt
) {

    try {
        if (opt.overwrite) { mpi.remove(opt.output); }
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

    HTS::_disable_logging();

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

    // Print the program title

    _print_title(mpi);
    __init_log(_LOG_FILE);

    // Delete the exiting output file if specified

    try {
        if (opt.overwrite) { mpi.remove(opt.output); }
    } catch (const std::exception& e) {
        mpi.err() << "Error: " << e.what() << "\n";
        return EXIT_FAILURE;
    }

    // Open the reference FASTA

    std::optional<BinaryFASTA> _fasta;
    try {

        if (mpi.root()) {
            _fasta.emplace(opt.fasta);
        }
        mpi.barrier();
        if (!mpi.root()) {
            _fasta.emplace(opt.fasta);
        }
        mpi.barrier();

    } catch (const std::exception& e) {
        mpi.err() << "Error: " << e.what() << "\n";
        __cleanup(mpi, opt);
        return EXIT_FAILURE;
    }
    BinaryFASTA fasta = std::move(_fasta.value());

    // Get the optimal chunksize given the number of processes

    int64_t _chunksize = mpi.chunksize(opt.chunk_size, fasta.size());

    // Create the HDF5 file

    std::optional<HDF5::File> _hdf5;
    try {
        _hdf5.emplace(opt.output, HDF5::RWC, mpi, _chunksize, opt.compression);
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

    size_t total = opt.files.value().size();
    if (total == 0) {
        if (opt.tokenize) {
            __write_sequences(fasta, hdf5, mpi);
        } else {
            mpi.out() << "        Nothing to do.\n";
        }
        mpi.down();
        return EXIT_SUCCESS;
    }

    // Open the HTS files, whose constructors are not thread-safe

    std::optional<HTS::FileGroup> _files;
    try {

        if (mpi.root()) {
            _files.emplace(opt.files);
        }
        mpi.barrier();
        if (!mpi.root()) {
            _files.emplace(opt.files);
        }
        mpi.barrier();

    } catch (const std::exception& e) {
        mpi.err() << "Error: " << e.what() << "\n";
        __cleanup(mpi, opt);
        return EXIT_FAILURE;
    }
    HTS::FileGroup files = std::move(_files.value());
    mpi.barrier();

    // Get the desired operation modes

    cmuts::Mode mode;
    cmuts::Spread spread;
    try {
        mode   = cmuts::mode(opt.lowmem, opt.joint);
        spread = cmuts::spread(opt.uniform_spread, opt.mutation_spread);
    } catch (const std::exception& e) {
        mpi.err() << "Error: " << e.what() << "\n";
        __cleanup(mpi, opt);
        return EXIT_FAILURE;
    }

    // Initialise the parameters for the main counting function

    cmuts::Params params = {
        mode,
        spread,
        opt.min_mapq,
        opt.min_quality,
        opt.min_length,
        opt.max_length,
        opt.max_indel_length,
        opt .quality_window,
        opt.collapse,
        !opt.no_mismatch,
        !opt.no_insertion,
        !opt.no_deletion,
        opt.subsample,
        opt.filter_coverage,
        !opt.disable_ambiguous,
        opt.contiguous_ambiguous,
        opt.print_every
    };

    // Initialise the stats tracker, and print the header

    cmuts::Stats stats(
        files.size(),
        files.aligned(),
        files.unaligned(),
        files.references(),
        fasta.longest(),
        mpi
    );
    stats.header();

    size_t processed = 0;

    for (auto& input : files) {

        std::unique_ptr<cmuts::Main> main;
        std::string name = _path(input->name());;

        try {
            main = cmuts::_get_main(*input, fasta, hdf5, mpi, params, stats, name);
            main->run();
            processed++;
        } catch (const std::exception& e) {
            mpi.err() << "Error processing the file \"" << input->name() << "\": " << e.what() << "\n";
            continue;
        }

    }

    if (processed == 0) {
        __cleanup(mpi, opt);
    } else if (opt.tokenize) {
        __write_sequences(fasta, hdf5, mpi);
    }

    mpi.down();

    if (processed < files.size()) {
        mpi.err() << "WARNING: only " << processed << " of " << files.size() << " files were processed.\n";
    }

    return EXIT_SUCCESS;

}
