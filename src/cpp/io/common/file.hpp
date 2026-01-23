#ifndef _CMUTS_IO_COMMON_FILE_HEADER
#define _CMUTS_IO_COMMON_FILE_HEADER

#include <cstdint>
#include <vector>
#include <string>
#include <memory>
#include <fstream>

extern "C" {
    #include <htslib/bgzf.h>
}

#include "alignment.hpp"


namespace HTS {


//
// Iterator
//


class Iterator {
protected:

    int64_t _reads = 0;
    int64_t _curr  = 0;

public:

    Iterator() = default;
    explicit Iterator(int64_t reads);
    virtual ~Iterator() = default;
    Iterator(const Iterator&) = delete;
    Iterator& operator=(const Iterator&) = delete;
    Iterator(Iterator&&) = delete;
    Iterator& operator=(Iterator&&) = delete;

    virtual Alignment next() = 0;
    bool end() const noexcept;

};


class EmptyIterator : public Iterator {
public:

    EmptyIterator() = default;
    Alignment next() override;

};


//
// Index
//


class IndexBlock {
public:

    int64_t ptr   = -1;
    int64_t reads = 0;

    bool empty() const;
    void tell(BGZF* _hts_bgzf);
    void write_ptr(std::ofstream& file);
    void write_bad_ptr(std::ofstream& file);
    void write_reads(std::ofstream& file);

};


class Index {
private:

    std::string _name;
    std::ifstream _file;
    int32_t _references;

public:

    Index() = default;
    Index(const std::string& filename, int32_t references);

    IndexBlock read(int32_t ix);
    int64_t aligned();
    int64_t unaligned();

};


Index _sam_index(BGZF* file);
Index _bam_index(BGZF* file);
Index _cram_index(BGZF* file);


//
// File
//


class File {
protected:

    std::string _name;
    FileType _type;

    BGZF* _hts_bgzf     = nullptr;
    bool _sorted        = false;
    int32_t _references = 0;

    Index _index;

public:

    explicit File(const std::string& name, FileType type);
    File(File&& other) noexcept;
    File& operator=(File&& other) = delete;
    File(const File&) = delete;
    File& operator=(const File&) = delete;
    virtual ~File();

    std::string name() const;
    FileType type() const;
    int32_t size() const;
    int64_t aligned();
    int64_t unaligned();

    virtual std::shared_ptr<Iterator> get(int32_t ix, bool seek) = 0;

};


std::unique_ptr<File> _get_sam(const std::string& name);
std::unique_ptr<File> _get_bam(const std::string& name);
std::unique_ptr<File> _get_cram(const std::string& name);


//
// FileGroup
//


class FileGroup {
private:

    std::vector<std::unique_ptr<File>> _group;

public:

    explicit FileGroup(const std::vector<std::string>& filenames);

    int32_t size() const;
    int32_t size(FileType type) const;
    int64_t aligned();
    int64_t unaligned();
    int32_t references() const;

    auto begin() -> decltype(_group.begin());
    auto begin() const -> decltype(_group.begin());
    auto end() -> decltype(_group.end());
    auto end() const -> decltype(_group.end());

};


} // namespace HTS


#endif
