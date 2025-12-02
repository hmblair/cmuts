#ifndef CMUTS_UTIL_HEADER
#define CMUTS_UTIL_HEADER

#include <stdexcept>
#include <string>
#include <iostream>
#include <fstream>
#include <ctime>
#include <execinfo.h>
#include <argparse/argparse.hpp>

#include "mpi.hpp"

const std::string _LOG_FILE = ".cmuts.log";
constexpr int32_t MAX_TRACE = 256;

typedef argparse::ArgumentParser Parser;

// Tag type to disambiguate group-only constructor from default value constructor
struct GroupTag {};

void __log(const std::string& filename, const std::string& message);
void __init_log(const std::string& filename);
[[noreturn]] void __throw_and_log(const std::string& filename, const std::string& err);
[[noreturn]] void __throw_and_log(
    const std::string& filename,
    const char* file,
    int line,
    const std::string& err
);

// Convenience macro for throwing errors with location information
#define CMUTS_THROW(msg) __throw_and_log(_LOG_FILE, __FILE__, __LINE__, msg)

// Convenience macro for conditional throws
#define CMUTS_THROW_IF(cond, msg) \
    do { if (cond) { CMUTS_THROW(msg); } } while (0)

//
// Structured logging with levels
//

enum class LogLevel { DEBUG = 0, INFO = 1, WARN = 2, ERROR = 3 };

// Global log level (default INFO, set via --verbose for DEBUG)
extern LogLevel g_log_level;

// Log function with level, file, and line context
void __log_level(
    const std::string& log_file,
    LogLevel level,
    const char* file,
    int line,
    const std::string& msg
);

// Set log level programmatically
void set_log_level(LogLevel level);
void set_log_level_from_verbose(bool verbose);

// Logging macros with file/line context
#define CMUTS_LOG(level, msg) \
    do { if (static_cast<int>(level) >= static_cast<int>(g_log_level)) \
        __log_level(_LOG_FILE, level, __FILE__, __LINE__, msg); } while(0)

#define CMUTS_DEBUG(msg) CMUTS_LOG(LogLevel::DEBUG, msg)
#define CMUTS_INFO(msg)  CMUTS_LOG(LogLevel::INFO, msg)
#define CMUTS_WARN(msg)  CMUTS_LOG(LogLevel::WARN, msg)
#define CMUTS_ERROR(msg) CMUTS_LOG(LogLevel::ERROR, msg)

//
// Validation macros
//

// Always-on validation (file I/O, initialization, user input)
#define CMUTS_CHECK_BOUNDS(idx, size, msg) \
    do { if (static_cast<size_t>(idx) >= static_cast<size_t>(size)) { \
        CMUTS_THROW(std::string(msg) + " (index=" + std::to_string(idx) + \
                    ", size=" + std::to_string(size) + ")"); \
    }} while(0)

#define CMUTS_CHECK_NOT_NULL(ptr, msg) \
    do { if ((ptr) == nullptr) { CMUTS_THROW(msg); }} while(0)

#define CMUTS_CHECK_POSITIVE(val, msg) \
    do { if ((val) <= 0) { CMUTS_THROW(std::string(msg) + " (value=" + std::to_string(val) + ")"); }} while(0)

// Debug-only validation (hot loops - disabled in release builds)
#ifdef NDEBUG
    #define CMUTS_DEBUG_CHECK_BOUNDS(idx, size, msg) ((void)0)
    #define CMUTS_DEBUG_CHECK_NOT_NULL(ptr, msg) ((void)0)
#else
    #define CMUTS_DEBUG_CHECK_BOUNDS(idx, size, msg) CMUTS_CHECK_BOUNDS(idx, size, msg)
    #define CMUTS_DEBUG_CHECK_NOT_NULL(ptr, msg) CMUTS_CHECK_NOT_NULL(ptr, msg)
#endif

template <typename T>
class Arg {
private:

    Parser& _parser;
    std::string _name;

public:

    Arg(
        Parser& parser,
        const std::string& short_name,
        const std::string& long_name,
        const std::string& help
    );
    Arg(
        Parser& parser,
        const std::string& short_name,
        const std::string& long_name,
        const std::string& help,
        GroupTag,
        const std::string& group
    );
    Arg(
        Parser& parser,
        const std::string& short_name,
        const std::string& long_name,
        const std::string& help,
        T default_value,
        const std::string& group = ""
    );
    T value() const;
    operator T() const;

};

class Program {
protected:

    Parser _parser;

public:

    Program(const std::string& program, const std::string& version);
    void parse(int argc, char** argv);

};

namespace Utils {

std::string _get_version();

void cursor_up(int num_lines);
void cursor_down(int num_lines);
void version();
void divider();

class Line {

private:

    std::string text;
    std::string suffix = "";
    int width          = 31;
    int precision      = 1;

public:

    Line(const std::string& text);
    Line(const std::string& text, const std::string& suffix);

    void print(const std::string& value) const;
    void print(int64_t value) const;
    void print(double value) const;

};

}

struct comma_out : std::numpunct<char>
{
    char do_thousands_sep()   const { return ','; }
    std::string do_grouping() const { return "\3"; }
};
void __imbue();

#endif
