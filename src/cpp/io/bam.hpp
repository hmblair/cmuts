#ifndef _CMUTS_BAM_HEADER
#define _CMUTS_BAM_HEADER

#include "common.hpp"
#include "infra/utils.hpp"

const char MD_DEL = '^';
const char MD_NULL = '0';

namespace HTS {

class MD_tag {
  private:
    // The MD tag

    std::string _tag;

    // The current position in the MD tag

    int32_t pos = 0;

    // The current op

    CIGAR_op curr;

    // Update curr to the op of the current position

    void update();

  public:
    // Create an MD_tag object from the raw MD tag

    explicit MD_tag(const std::string& tag);

    // Return the next operation in the MD tag, possibly truncated to a length of max.
    // Store the remainder of the truncation in curr.

    CIGAR_op advance();
    CIGAR_op advance(int32_t max);

    // The MD tag as it would appear in a SAM file

    std::string str() const;
};

//
// bamAlignmentIter
//

class BamIterator : public Iterator {
  private:
    BGZF* _hts_bgzf = nullptr;
    bam1_t* _hts_aln = nullptr;

  public:
    BamIterator(BGZF* _hts_bgzf, bam1_t* _hts_aln, int64_t reads);
    BamIterator(BamIterator&& other) noexcept;
    BamIterator& operator=(BamIterator&& other) noexcept;
    BamIterator(BamIterator& other) = delete;
    BamIterator operator=(BamIterator& other) = delete;
    ~BamIterator() override = default;

    Alignment next() override;
};

//
// BamFile
//

class BamFile : public File {
  private:
    bam1_t* _hts_aln = nullptr;

  public:
    explicit BamFile(const std::string& name);
    ~BamFile() override;

    std::shared_ptr<Iterator> get(int32_t ix, bool seek) override;
};

//
// SamIterator
//
// SAM records decode to the same bam1_t as BAM, so the only difference from
// BamIterator is how the next record is obtained: a line of text is read from
// the (transparently uncompressed) BGZF stream and parsed with sam_parse1.
//

class SamIterator : public Iterator {
  private:
    BGZF* _hts_bgzf = nullptr;
    bam1_t* _hts_aln = nullptr;
    sam_hdr_t* _hdr = nullptr;
    kstring_t _line = KS_INITIALIZE;

  public:
    SamIterator(BGZF* _hts_bgzf, bam1_t* _hts_aln, sam_hdr_t* _hdr, int64_t reads);
    SamIterator(SamIterator&&) = delete;
    SamIterator& operator=(SamIterator&&) = delete;
    SamIterator(const SamIterator&) = delete;
    SamIterator& operator=(const SamIterator&) = delete;
    ~SamIterator() override;

    Alignment next() override;
};

//
// SamFile
//

class SamFile : public File {
  private:
    bam1_t* _hts_aln = nullptr;
    sam_hdr_t* _hdr = nullptr;

  public:
    explicit SamFile(const std::string& name);
    ~SamFile() override;

    std::shared_ptr<Iterator> get(int32_t ix, bool seek) override;
};

} // namespace HTS

#endif
