#ifndef _CMUTS_OPTIONS_HEADER_
#define _CMUTS_OPTIONS_HEADER_

#include <string>
#include <vector>

namespace cmuts {

struct OptionConfig {
    std::string short_name;
    std::string long_name;
    std::string help;
};

template<typename T>
struct TypedOptionConfig : OptionConfig {
    T default_value;
    TypedOptionConfig(const std::string& s, const std::string& l, const std::string& h, const T& def = T{})
        : OptionConfig{s, l, h}, default_value(def) {}
};

namespace options {

    // Input/Output options
    inline const OptionConfig FILES{"", "files", "The input SAM/BAM/CRAM files."};
    inline const OptionConfig OUTPUT{"-o", "--output", "The output HDF5 file."};
    inline const OptionConfig FASTA{"-f", "--fasta", "The reference FASTA file."};
    inline const OptionConfig OVERWRITE{"", "--overwrite", "Overwrite an existing HDF5 file."};
    inline const OptionConfig REBUILD{"", "--rebuild", "Rebuild all index files."};

    // Quality and filtering options
    inline const TypedOptionConfig<int> MIN_PHRED{"", "--min-phred", "PHRED score threshold for base processing.", 10};
    inline const TypedOptionConfig<int> MIN_MAPQ{"", "--min-mapq", "Mapping quality threshold for alignment processing.", 10};
    inline const TypedOptionConfig<int> MIN_LENGTH{"", "--min-length", "Skip reads shorter than this length.", 2};
    inline const TypedOptionConfig<int> MAX_LENGTH{"", "--max-length", "Skip reads longer than this length.", 1024};
    inline const TypedOptionConfig<int> MAX_HAMMING{"", "--max-hamming", "The maximum number of mismatches, insertions, and deletions in a processed read.", 1024};
    inline const TypedOptionConfig<float> SUBSAMPLE{"", "--subsample", "Randomly choose to use a read with this probability.", 1.0f};

    // Processing options
    inline const TypedOptionConfig<int> COMPRESSION{"-c", "--compression", "Compression level of the HDF5 output (0-9).", 3};
    inline const TypedOptionConfig<int> MAX_INDEL_LENGTH{"", "--max-indel-length", "The longest indels to consider.", 10};
    inline const TypedOptionConfig<int> CHUNK_SIZE{"", "--chunk-size", "The number of references to process at a time per thread.", 128};
    inline const TypedOptionConfig<int> QUALITY_WINDOW{"", "--quality-window", "Check the quality of each base in a window of this size around each base.", 2};
    inline const TypedOptionConfig<int> COLLAPSE{"", "--collapse", "Collapse modifications within this distance of each other in a given read.", 2};

    // Mode options
    inline const OptionConfig JOINT{"", "--pairwise", "Compute pairwise modification counts."};
    inline const OptionConfig TOKENIZE{"", "--tokenize", "Tokenize the reference sequences."};

    // Mutation type filters
    inline const OptionConfig NO_MISMATCH{"", "--no-mismatches", "Do not count mismatches as modifications."};
    inline const OptionConfig NO_INSERTION{"", "--no-insertions", "Do not count insertions as modifications."};
    inline const OptionConfig NO_DELETION{"", "--no-deletions", "Do not count deletions as modifications."};

    // Strand options
    inline const OptionConfig NO_REVERSE{"", "--no-reverse", "Ignore reverse-complemented reads."};
    inline const OptionConfig ONLY_REVERSE{"", "--only-reverse", "Use only reverse-complemented reads."};

    // Deletion handling
    inline const OptionConfig UNIFORM_SPREAD{"", "--uniform-spread", "Uniformly spread out ambiguous deletions."};
    inline const OptionConfig NO_SPREAD{"", "--no-spread", "Do not spread ambiguous deletions."};
    inline const TypedOptionConfig<int> DELETION_GAP{"", "--deletion-gap", "The number of gaps to allow when detecting ambiguous deletions.", 0};
    inline const OptionConfig DISABLE_AMBIGUOUS{"", "--disable-ambiguous", "Disable the ambiguous delection detection algorithm, relying on the deletion provided by the alignment."};

    // Quality filtering
    inline const OptionConfig NO_FILTER_MATCHES{"", "--no-match-filter", "Do not filter matches based on their PHRED base score."};
    inline const OptionConfig NO_FILTER_INSERTIONS{"", "--no-insertion-filter", "Do not filter insertions based on their PHRED base score."};
    inline const OptionConfig NO_FILTER_DELETIONS{"", "--no-deletion-filter", "Do not filter deletions based on their PHRED base score."};

    // Base filtering
    inline const TypedOptionConfig<std::string> IGNORE_BASES{"", "--ignore-bases", "Do not count mismatches or deletions occuring at these bases. Pass as a single string", ""};

}

} // namespace cmuts

#endif
