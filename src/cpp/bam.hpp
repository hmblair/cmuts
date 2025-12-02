# ifndef _CMUTS_BAM_HEADER
# define _CMUTS_BAM_HEADER

#include "common.hpp"
#include "utils.hpp"

const char MD_DEL  = '^';
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

    BGZF* _hts_bgzf  = nullptr;
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





} // namespace HTS





#endif
