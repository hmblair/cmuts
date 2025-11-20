#include "bam.hpp"
#include "common.hpp"
#include "htslib/hts.h"
#include <stdexcept>





namespace HTS {





//
// Converting to and from internal HTSlib representaitons
//





static inline base_t _base_from_hts(uint8_t base) {

    switch (base) {
        case HTS_A: { return IX_A;   }
        case HTS_C: { return IX_C;   }
        case HTS_G: { return IX_G;   }
        case HTS_T: { return IX_T;   }
        default:    { return IX_UNK; }
    }

}


static inline CIGAR_op _cigar_from_hts(uint32_t op) {

    int32_t length = bam_cigar_oplen(op);
    CIGAR_t type   = CIGAR_t::UNKNOWN;

    switch (bam_cigar_op(op)) {
        case BAM_CMATCH:
        case BAM_CEQUAL: {
            type = CIGAR_t::MATCH;
            break;
        }
        case BAM_CDIFF: {
            type = CIGAR_t::MISMATCH;
            break;
        }
        case BAM_CDEL: {
            type = CIGAR_t::DEL;
            break;
        }
        case BAM_CINS: {
            type = CIGAR_t::INS;
            break;
        }
        case BAM_CSOFT_CLIP: {
            type = CIGAR_t::SOFT;
            break;
        }
        case BAM_CHARD_CLIP: {
            type = CIGAR_t::HARD;
            break;
        }
        case BAM_CREF_SKIP: {
            type = CIGAR_t::SKIP;
            break;
        }
        case BAM_CPAD: {
            type = CIGAR_t::PAD;
            break;
        }
        default: {
            type = CIGAR_t::UNKNOWN;
            break;
        }
    }

    return {type, length};

}


static inline base_t _bam_base(bam1_t* _hts_aln, int32_t ix) {

    uint8_t* _query   = bam_get_seq(_hts_aln);
    uint8_t _hts_base = bam_seqi(_query, ix);

    return _base_from_hts(_hts_base);

}





//
// MD tag
//





static inline bool _is_digit(const std::string& tag, int32_t pos) {

    return (
	pos < static_cast<int32_t>(tag.size()) &&
	static_cast<bool>(std::isdigit(tag[pos]))
    );

}


static inline bool _is_alpha(const std::string& tag, int32_t pos) {

    return (
        pos < static_cast<int32_t>(tag.size()) &&
        static_cast<bool>(std::isalpha(tag[pos]))
    );

}


static inline bool _is_deletion(const std::string& tag, int32_t pos) {

    return pos < static_cast<int32_t>(tag.size()) && tag[pos] == MD_DEL;

}


static inline bool _is_null(const std::string& tag, int32_t pos) {

    return pos < static_cast<int32_t>(tag.size()) && tag[pos] == MD_NULL;

}


static inline void _skip_null(const std::string& tag, int32_t& pos) {

    if (_is_null(tag, pos)) { pos++; }

}


static inline int32_t _count_matches(const std::string& tag, int32_t& pos) {

    static const int32_t BASE = 10;
    int32_t _matches = 0;

    while (_is_digit(tag, pos)) {
        _matches = _matches * BASE + (tag[pos] - '0');
        pos++;
    }

    return _matches;

}


static inline int32_t _count_mismatches(const std::string& tag, int32_t& pos, base_t& base) {

    int32_t _mismatches = 0;
    while (_is_alpha(tag, pos)) {
        base = _from_char(tag[pos]);
        _mismatches++;
        pos++;
        _skip_null(tag, pos);
    }

    return _mismatches;

}


static inline int32_t _count_deletions(const std::string& tag, int32_t& pos) {

    int32_t _deletions = 0;
    pos++;
    while (_is_alpha(tag, pos)) {
        _deletions++;
        pos++;
    }

    _skip_null(tag, pos);
    return _deletions;

}


static inline CIGAR_t _md_type(const std::string& tag, int32_t pos) {

    if (_is_digit(tag, pos) && !_is_null(tag, pos)) {
        return CIGAR_t::MATCH;
    } else if (_is_alpha(tag, pos)) {
        return CIGAR_t::MISMATCH;
    } else if (_is_deletion(tag, pos)) {
        return CIGAR_t::DEL;
    } else {
        return CIGAR_t::UNKNOWN;
    }

}


MD_tag::MD_tag(const std::string& tag) : _tag(tag) {

    _skip_null(_tag, pos);

}


std::string MD_tag::str() const {

    return _tag;

}


void MD_tag::update() {

    CIGAR_t type  = _md_type(_tag, pos);
    int32_t count = 0;
    base_t base   = IX_UNK;

    switch (type) {

        case CIGAR_t::MATCH: {
            count = _count_matches(_tag, pos);
            break;
        }

        case CIGAR_t::MISMATCH: {
            count = _count_mismatches(_tag, pos, base);
            break;
        }

        case CIGAR_t::DEL: {
            count = _count_deletions(_tag, pos);
            break;
        }

        default: {
            __throw_and_log(_LOG_FILE, "Unknown MD tag operation \"" + std::to_string(_tag[pos]) + "\" at tag position " + std::to_string(pos) + ". The tag is " + _tag + ".");
        }

    }

    curr = CIGAR_op(type, count);
    curr.base = base;

}


CIGAR_op MD_tag::advance() {

    // If the current op is empty, then move to the next op
    if (curr.empty()) { update(); }

    // Return the current op and set it to be empty
    CIGAR_op _tmp = curr;
    curr = CIGAR_op();
    return _tmp;

}


CIGAR_op MD_tag::advance(int32_t max) {

    // If the current op is empty, then move to the next op
    if (curr.empty()) { update(); }

    // Split the op based on the requested maximum size;
    // return the first half and store the second half
    std::pair<CIGAR_op, CIGAR_op> pair = curr.split(max);
    curr = pair.second;
    return pair.first;

}


static inline CIGAR _get_cigar_core(
    std::span<const uint32_t> hts_cigar,
    const std::string& md_tag_str,
    bam1_t* _hts_aln
) {

    int32_t qpos = 0;

    MD_tag md_tag(md_tag_str);

    // Reserve the minimum amount of space required

    auto size = static_cast<int32_t>(hts_cigar.size() + 1);
    CIGAR cigar(size);

    // Add the termination event first

    CIGAR_op term(CIGAR_t::TERM);
    cigar.append(term);

    // Read the remaining operations from the HTS CIGAR string

    for (const auto& hts_op : hts_cigar) {

        CIGAR_op op = _cigar_from_hts(hts_op);

        // If the CIGAR string indicates a match, we must use the MD tag
        // in order to find any mismatches

        if (op.type() == CIGAR_t::MATCH) {

            int32_t remaining = op.length();

            while (remaining > 0) {

                CIGAR_op _md_cig = md_tag.advance(remaining);
                qpos += (_md_cig.qlength() - 1);

                if (_md_cig.type() == CIGAR_t::MISMATCH) {
                    _md_cig.base = _bam_base(_hts_aln, qpos);
                }

                cigar.append(_md_cig);
                remaining -= _md_cig.length();
                qpos++;

            }

        }

        // Else, the HTS CIGAR string contains all the information we need.
        // The MD tag must be advanced if it is a deletion to keep the two in sync.

        else {

            qpos += (op.qlength() - 1);

            if (op.type() == CIGAR_t::INS) {
                op.base = _bam_base(_hts_aln, qpos);
            }

            if (op.type() == CIGAR_t::DEL) {
                (void)md_tag.advance();
            }

            cigar.append(op);
            qpos++;

        }


    }

    return cigar;

}





//
// BAM-related functions
//





static inline bam1_t* _open_aln() {

    bam1_t* _hts_aln = bam_init1();
    if (_hts_aln == nullptr) {
        __throw_and_log(_LOG_FILE, "Failed to allocate memory for a sequence aligment.");
    }

    // Set the empty alignment to be unaligned

    _hts_aln->core.tid = -1;
    return _hts_aln;

}


static inline void _close_aln(bam1_t*& _hts_aln) {

    if (_hts_aln != nullptr) {
        bam_destroy1(_hts_aln);
        _hts_aln = nullptr;
    }

}


static inline void _read_bam(BGZF* _hts_bgzf, bam1_t* _hts_aln) {

    if (bam_read1(_hts_bgzf, _hts_aln) < 0) {
        __throw_and_log(_LOG_FILE, "The alignment ended prematurely.");
    }

}


static inline bool _bam_aligned(bam1_t* _hts_aln) {

    return _hts_aln->core.tid >= 0;

}


static inline bool _bam_primary(bam1_t* _hts_aln) {

    int32_t mask = 0x100;
    return !static_cast<bool>(_hts_aln->core.flag & mask);

}


static inline int32_t _bam_length(bam1_t* _hts_aln) {

    return _hts_aln->core.l_qseq;

}


static inline int32_t _bam_offset(bam1_t* _hts_aln) {

    return _hts_aln->core.pos;

}


static inline qual_t _bam_mapq(bam1_t* _hts_aln) {

    return _hts_aln->core.qual;

}


static inline std::string _bam_md(bam1_t* _hts_aln) {

    uint8_t* md_ptr = bam_aux_get(_hts_aln, "MD");
    if (md_ptr == nullptr) {
        __throw_and_log(_LOG_FILE, "Failed to retrieve the MD tag.");
    }
    char* tag = bam_aux2Z(md_ptr);

    return {tag, strlen(tag)};

}


static inline std::span<const uint32_t> _bam_cigar_str(bam1_t* _hts_aln) {

    uint32_t* raw_str = bam_get_cigar(_hts_aln);
    if (raw_str == nullptr) {
        __throw_and_log(_LOG_FILE, "Error retrieving the CIGAR string.");
    }
    hts_pos_t length = _hts_aln->core.n_cigar;

    return {raw_str, raw_str + length};

}


static inline CIGAR _bam_cigar(bam1_t* _hts_aln) {
    std::string _md_tag = _bam_md(_hts_aln);
    std::span<const uint32_t> _cigar_str = _bam_cigar_str(_hts_aln);
    return _get_cigar_core(_cigar_str, _md_tag, _hts_aln);

}


static inline PHRED _bam_phred(bam1_t* _hts_aln) {

    qual_t* _quality = bam_get_qual(_hts_aln);
    int32_t length = _bam_length(_hts_aln);
    std::vector<qual_t> scores = {_quality, _quality + length};
    return PHRED(scores);

}


bamIterator::bamIterator(BGZF* _hts_bgzf, bam1_t* _hts_aln, int64_t reads)
    : Iterator(reads), _hts_bgzf(_hts_bgzf), _hts_aln(_hts_aln) {}


bamIterator::bamIterator(bamIterator&& other) noexcept
    : Iterator(other._reads),
      _hts_bgzf(other._hts_bgzf),
      _hts_aln(other._hts_aln) {

    other._hts_bgzf = nullptr;
    other._hts_aln  = nullptr;

}


bamIterator& bamIterator::operator=(bamIterator&& other) noexcept {

    if (this != &other) {
        _hts_bgzf = other._hts_bgzf;
        _hts_aln  = other._hts_aln;
        _reads = other._reads;
        _curr = other._curr;
        other._hts_bgzf = nullptr;
        other._hts_aln  = nullptr;
    }

    return *this;

}


Alignment bamIterator::next() {

    _read_bam(_hts_bgzf, _hts_aln);
    _curr++;

    Alignment aln;

    aln.aligned  = _bam_aligned(_hts_aln);
    aln.primary  = _bam_primary(_hts_aln);
    aln.reversed = bam_is_rev(_hts_aln);
    aln.mapq     = _bam_mapq(_hts_aln);
    aln.length   = _bam_length(_hts_aln);
    aln.offset   = _bam_offset(_hts_aln);
    aln.cigar    = _bam_cigar(_hts_aln);
    aln.phred    = _bam_phred(_hts_aln);

    return aln;

}





//
// bamFile
//





static inline Header _read_bam_header(BGZF* _bgzf_file) {

    FileType type = _get_filetype(_bgzf_file);
    if (type != FileType::BAM) {
        __throw_and_log(_LOG_FILE, "Opening a non-BAM file using _read_bam_header.");
    }

    auto length = _read_bgzf_single<int32_t>(_bgzf_file);
    std::unique_ptr<ByteStream> stream = std::make_unique<bgzfFileStream>(_bgzf_file, length);

    // Skip the initial, uncompressed SAM header

    Header data = _read_sam_header(stream, false);

    // Read the number of references

    _read_bgzf(_bgzf_file, &data.references, sizeof(uint32_t));

    // Skip the remainder of the compressed portion

    char buffer[BGZF_BUFFER + sizeof(uint32_t)];
    uint32_t name_len = 0;

    for (int32_t ix = 0; ix < data.references; ix++) {
        _read_bgzf(_bgzf_file, &name_len, sizeof(uint32_t));
        if (name_len >= BGZF_BUFFER) {
            __throw_and_log(_LOG_FILE, "Header name larger than the buffer.");
        }
        _read_bgzf(_bgzf_file, buffer, name_len + sizeof(uint32_t));
    }

    return data;

}


static inline void _build_bam_index(
    BGZF* _bgzf_file,
    const std::string& filename,
    int32_t references
) {

    _throw_if_exists(filename);
    std::ofstream outfile(filename);

    bam1_t* _hts_aln = _open_aln();
    IndexBlock block;
    int64_t unaligned = 0;

    // The beginning of the non-header section

    block.tell(_bgzf_file);
    outfile.seekp(0);
    block.write_ptr(outfile);

    int32_t tid = 0, curr;
    while (bam_read1(_bgzf_file, _hts_aln) > 0) {

        curr = _hts_aln->core.tid;
        if (curr > tid) {

            block.write_reads(outfile);
            block.reads = 0;

            // Write empty blocks for all references which did not appear in the file

            for (int32_t jx = tid + 1; jx < curr; jx++) {
                block.write_ptr(outfile);
                block.write_reads(outfile);
            }

            block.reads++;
            tid = curr;
            block.write_ptr(outfile);

        } else if (curr == -1) {

            unaligned++;

        } else {

            block.reads++;

        }

        block.tell(_bgzf_file);

    }

    block.write_reads(outfile);
    block.reads = 0;

    // Write any final references with no aligned reads

    for (int32_t jx = tid + 1; jx < references; jx++) {
        block.write_ptr(outfile);
        block.write_reads(outfile);
    }

    // Write the unaligned read count

    outfile.write(reinterpret_cast<char*>(&unaligned), sizeof(int64_t));
    outfile.close();

}


bamFile::bamFile(const std::string& name)
    : File(name, FileType::BAM), _hts_aln(_open_aln()) {

    Header header = _read_bam_header(_hts_bgzf);
    _references = header.references;

    std::string _index_name = _name + CMUTS_INDEX;
    if (!std::filesystem::exists(_index_name) && !cmuts::mutex::check(_name)) {
        _build_bam_index(_hts_bgzf, _index_name, _references);
        __log(_LOG_FILE, "Successfully created " + _index_name + ".");
    }

    // Will trigger if another thread started creating the index file first

    cmuts::mutex::wait(_name);

    _index = Index(_index_name, _references);
    __log(_LOG_FILE, "Successfully loaded " + _index_name + ".");

}


bamFile::~bamFile() {

    _close_aln(_hts_aln);

}


std::shared_ptr<Iterator> bamFile::get(int32_t ix, bool seek) {

    IndexBlock block = _index.read(ix);
    if (seek) { _seek_bgzf(_hts_bgzf, block.ptr); }
    return std::make_shared<bamIterator>(_hts_bgzf, _hts_aln, block.reads);

}





//
// From common.hpp
//





std::unique_ptr<File> _get_sam(const std::string& name) {

    return std::make_unique<bamFile>(name);

}


std::unique_ptr<File> _get_bam(const std::string& name) {

    return std::make_unique<bamFile>(name);

}




} // namespace HTS
