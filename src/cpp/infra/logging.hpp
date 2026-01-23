#ifndef CMUTS_LOGGING_HEADER
#define CMUTS_LOGGING_HEADER

#include <stdexcept>
#include <string>

// Log file name
const std::string _LOG_FILE = ".cmuts.log";

// Maximum stack trace depth
constexpr int32_t MAX_TRACE = 256;

// Basic logging functions
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
// Configure via environment variables:
//   CMUTS_LOG_LEVEL: TRACE, DEBUG, INFO, WARN, ERROR (default: INFO)
//   CMUTS_LOG_STDERR: 1 to enable real-time stderr output (default: off)
//

enum class LogLevel { TRACE = 0, DEBUG = 1, INFO = 2, WARN = 3, ERROR = 4 };

// Global log level (default INFO, configurable via env or --verbose)
extern LogLevel g_log_level;

// Whether to also output to stderr (for real-time debugging)
extern bool g_log_stderr;

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

// Initialize logging from environment variables
// Call early in main(), before other logging calls
void init_logging_from_env();

// Logging macros with file/line context
#define CMUTS_LOG(level, msg) \
    do { if (static_cast<int>(level) >= static_cast<int>(g_log_level)) \
        __log_level(_LOG_FILE, level, __FILE__, __LINE__, msg); } while(0)

#define CMUTS_TRACE(msg) CMUTS_LOG(LogLevel::TRACE, msg)
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

#endif
