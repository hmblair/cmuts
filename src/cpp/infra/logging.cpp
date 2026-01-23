#include "infra/logging.hpp"

#include <chrono>
#include <cstdlib>
#include <ctime>
#include <execinfo.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

#ifdef MPI_BUILD
#include <mpi.h>
#endif

#ifndef VERSION
#define VERSION "unknown"
#endif

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
