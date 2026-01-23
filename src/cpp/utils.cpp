#include "utils.hpp"

#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <sstream>

static inline void __log_trace(std::ofstream& file) {

    void* array[MAX_TRACE];
    int size = backtrace(array, MAX_TRACE);
    char** symbols = backtrace_symbols(array, size);

    for (size_t i = 0; i < size; i++) {
        file << symbols[i] << "\n";
    }

    free(symbols);

}

void __log(const std::string& filename, const std::string& message) {

    std::ofstream file(filename, std::ios::app);

    if (!file) {
        std::cerr << "Error opening the log file \"" << filename << "\".\n";
        return;
    }

    file << "MESSAGE: " << message;

    #ifdef MPI_BUILD
        int rank = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        file << " (process " << std::to_string(rank) << ")";
    #endif

    file << std::endl;

}

void __init_log(const std::string& filename) {

    std::time_t now = std::time(nullptr);
    std::ofstream file(filename, std::ios::app);

    if (!file) {
        std::cerr << "Error opening the log file \"" << filename << "\".\n";
        return;
    }

    #ifdef MPI_BUILD
        int rank = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank == 0) {
    #endif

    file << "--------------- " << ("cmuts v" + std::string(VERSION) + " - ") << std::ctime(&now);

    #ifdef MPI_BUILD
        }
    #endif

}

[[noreturn]] void __throw_and_log(const std::string& filename, const std::string& err) {

    std::ofstream file(filename, std::ios::app);

    if (!file) {
        std::cerr << "Error opening the log file \"" << filename << "\".\n";
    } else {
        file << "ERROR: " << err << std::endl;
        __log_trace(file);
    }

    throw std::runtime_error(err);

}

[[noreturn]] void __throw_and_log(
    const std::string& filename,
    const char* source_file,
    int line,
    const std::string& err
) {
    std::ofstream file(filename, std::ios::app);

    if (!file) {
        std::cerr << "Error opening the log file \"" << filename << "\".\n";
    } else {
        file << "ERROR [" << source_file << ":" << line << "]: " << err << std::endl;
        __log_trace(file);
    }

    // Include location in the exception message for easier debugging
    std::string full_msg = err + " [" + source_file + ":" + std::to_string(line) + "]";
    throw std::runtime_error(full_msg);
}

//
// Structured logging
//

LogLevel g_log_level = LogLevel::INFO;
bool g_log_stderr = false;

static const char* _level_names[] = {"TRACE", "DEBUG", "INFO", "WARN", "ERROR"};

// Get elapsed time since program start (thread-safe)
static double __get_elapsed_seconds() {
    static auto start_time = std::chrono::steady_clock::now();
    auto now = std::chrono::steady_clock::now();
    return std::chrono::duration<double>(now - start_time).count();
}

// Format a log line with all context
static std::string __format_log_line(
    LogLevel level,
    const char* file,
    int line,
    const std::string& msg
) {
    std::ostringstream oss;

    // Timestamp: elapsed seconds since start
    double elapsed = __get_elapsed_seconds();
    oss << std::fixed << std::setprecision(3) << "[" << elapsed << "s] ";

    // Level
    oss << "[" << _level_names[static_cast<int>(level)] << "] ";

    // MPI rank
    #ifdef MPI_BUILD
        int rank = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        oss << "[P" << rank << "] ";
    #endif

    // Source location (basename only for readability)
    const char* basename = file;
    for (const char* p = file; *p; ++p) {
        if (*p == '/') basename = p + 1;
    }
    oss << basename << ":" << line << " - " << msg;

    return oss.str();
}

void __log_level(
    const std::string& log_file,
    LogLevel level,
    const char* file,
    int line,
    const std::string& msg
) {
    std::string formatted = __format_log_line(level, file, line, msg);

    // Write to file
    std::ofstream log(log_file, std::ios::app);
    if (log) {
        log << formatted << "\n";
    }

    // Write to stderr if enabled (for real-time debugging)
    if (g_log_stderr) {
        std::cerr << formatted << std::endl;
    }
}

void set_log_level(LogLevel level) {
    g_log_level = level;
}

void set_log_level_from_verbose(bool verbose) {
    // Only upgrade to DEBUG if not already set lower by env var
    if (verbose && g_log_level > LogLevel::DEBUG) {
        g_log_level = LogLevel::DEBUG;
    }
}

void init_logging_from_env() {
    // Parse CMUTS_LOG_LEVEL
    const char* level_env = std::getenv("CMUTS_LOG_LEVEL");
    if (level_env) {
        std::string level_str(level_env);
        // Convert to uppercase for comparison
        for (auto& c : level_str) c = std::toupper(c);

        if (level_str == "TRACE") {
            g_log_level = LogLevel::TRACE;
        } else if (level_str == "DEBUG") {
            g_log_level = LogLevel::DEBUG;
        } else if (level_str == "INFO") {
            g_log_level = LogLevel::INFO;
        } else if (level_str == "WARN" || level_str == "WARNING") {
            g_log_level = LogLevel::WARN;
        } else if (level_str == "ERROR") {
            g_log_level = LogLevel::ERROR;
        }
    }

    // Parse CMUTS_LOG_STDERR
    const char* stderr_env = std::getenv("CMUTS_LOG_STDERR");
    if (stderr_env && (std::string(stderr_env) == "1" || std::string(stderr_env) == "true")) {
        g_log_stderr = true;
    }
}

template <typename T>
static inline void _add_arg(
    Parser& parser,
    const std::string& short_name,
    const std::string& long_name,
    const std::string& help
);

template <typename T>
static inline void _add_arg(
    Parser& parser,
    const std::string& short_name,
    const std::string& long_name,
    const std::string& help,
    T default_value
);

template<>
void _add_arg<bool>(
    Parser& parser,
    const std::string& short_name,
    const std::string& long_name,
    const std::string& help
) {
    parser.add_argument(short_name, long_name)
        .default_value(false)
        .implicit_value(true)
        .help(help);
}

template<>
void _add_arg<bool>(
    Parser& parser,
    const std::string& short_name,
    const std::string& long_name,
    const std::string& help,
    bool default_value
) {
    parser.add_argument(short_name, long_name)
        .default_value(default_value)
        .implicit_value(!default_value)
        .help(help);
}

template<>
void _add_arg<int>(
    Parser& parser,
    const std::string& short_name,
    const std::string& long_name,
    const std::string& help
) {
    parser.add_argument(short_name, long_name)
        .required()
        .scan<'d', int>()
        .help(help);
}

template<>
void _add_arg<int>(
    Parser& parser,
    const std::string& short_name,
    const std::string& long_name,
    const std::string& help,
    int default_value
) {
    parser.add_argument(short_name, long_name)
        .scan<'d', int>()
        .default_value(default_value)
        .help(help);
}

template<>
void _add_arg<float>(
    Parser& parser,
    const std::string& short_name,
    const std::string& long_name,
    const std::string& help
) {
    parser.add_argument(short_name, long_name)
        .required()
        .scan<'g', float>()
        .help(help);
}

template<>
void _add_arg<float>(
    Parser& parser,
    const std::string& short_name,
    const std::string& long_name,
    const std::string& help,
    float default_value
) {
    parser.add_argument(short_name, long_name)
        .scan<'g', float>()
        .default_value(default_value)
        .help(help);
}

template<>
void _add_arg<std::string>(
    Parser& parser,
    const std::string& short_name,
    const std::string& long_name,
    const std::string& help
) {
    parser.add_argument(short_name, long_name)
        .required()
        .help(help);
}

template<>
void _add_arg<std::string>(
    Parser& parser,
    const std::string& short_name,
    const std::string& long_name,
    const std::string& help,
    std::string default_value
) {
    parser.add_argument(short_name, long_name)
        .default_value(default_value)
        .help(help);
}

template<>
void _add_arg<std::vector<std::string>>(
    Parser& parser,
    const std::string& short_name,
    const std::string& long_name,
    const std::string& help
) {
    parser.add_argument(short_name, long_name)
        .remaining()
        .required()
        .help(help);
}

template<>
void _add_arg<std::vector<std::string>>(
    Parser& parser,
    const std::string& short_name,
    const std::string& long_name,
    const std::string& help,
    std::vector<std::string> default_value
) {
    parser.add_argument(short_name, long_name)
        .remaining()
        .default_value(default_value)
        .help(help);
}

template <typename T>
Arg<T>::Arg (
    Parser& parser,
    const std::string& short_name,
    const std::string& long_name,
    const std::string& help
) : _parser(parser), _name(long_name) {
    _add_arg<T>(parser, short_name, long_name, help);
}

template <typename T>
Arg<T>::Arg (
    Parser& parser,
    const std::string& short_name,
    const std::string& long_name,
    const std::string& help,
    GroupTag,
    const std::string& group
) : _parser(parser), _name(long_name) {
    if (!group.empty()) { parser.add_group(group); }
    _add_arg<T>(parser, short_name, long_name, help);
}

template <typename T>
Arg<T>::Arg (
    Parser& parser,
    const std::string& short_name,
    const std::string& long_name,
    const std::string& help,
    T default_value,
    const std::string& group
) : _parser(parser), _name(long_name) {
    if (!group.empty()) { parser.add_group(group); }
    _add_arg<T>(parser, short_name, long_name, help, default_value);
}

template <typename T>
T Arg<T>::value() const {
    return _parser.template get<T>(_name);
}

template <typename T>
Arg<T>::operator T() const {
    return value();
}

template class Arg<int>;
template class Arg<float>;
template class Arg<bool>;
template class Arg<std::string>;
template class Arg<std::vector<std::string>>;

Program::Program(const std::string& program, const std::string& version)
    : _parser(program, version) {}

void Program::parse(int argc, char** argv) {
    _parser.parse_args(argc, argv);
}

namespace Utils {

void cursor_up(int num_lines) {
    std::cout << "\033[" + std::to_string(num_lines) + "A";
}

void cursor_down(int num_lines) {
    for (int ix = 0; ix < num_lines; ix++) {
        std::cout << "\n";
    }
}

std::string SPACE = " ";
std::string DIV   = "─";

static inline std::string __repeat(const std::string& str, size_t count) {

    std::string result;
    result.reserve(str.size() * count);
    for (size_t ix = 0; ix < count; ++ix) { result += str; }
    return result;

}

#ifdef MPI_BUILD
#ifdef DEBUG
const std::string PROGRAM = "cmuts core MPI (DEBUG)";
#else
const std::string PROGRAM = "cmuts core MPI";
#endif
#else
#ifdef DEBUG
const std::string PROGRAM = "cmuts core (DEBUG)";
#else
const std::string PROGRAM = "cmuts core";
#endif
#endif

void version() {
    std::cout << __repeat(SPACE, 8) + PROGRAM + " version " + VERSION + "\n";
}

void divider() {
    std::cout << __repeat(SPACE, 6) + __repeat(DIV, 35) + "\n";
}

Line::Line(const std::string& text) : text(text) {}

Line::Line(const std::string& text, const std::string& suffix) : text(text), suffix(suffix) {}

void Line::print(const std::string& value) const {
    int spacing = width - text.length() - suffix.length();
    std::cout << "        " << text << ":"
              << std::setw(spacing) << value
              << suffix << "\n";
}

void Line::print(int64_t value) const {
    int spacing = width - text.length() - suffix.length();
    std::cout << "        " << text << ":"
              << std::setw(spacing) << value
              << suffix << "\n";
}

void Line::print(double value) const {
    int spacing = width - text.length() - suffix.length();
    std::cout << std::fixed << std::setprecision(precision)
              << "        " << text << ":"
              << std::setw(spacing) << value
              << suffix << "\n";
}

} // namepsace Utils

void __imbue() {
    std::cout.imbue(std::locale(std::cout.getloc(), new comma_out));
};
