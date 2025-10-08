#include "utils.hpp"
#include <stdexcept>

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

    file << "--------------- " << std::ctime(&now);

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
    T default_value
) : _parser(parser), _name(long_name) {
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

std::filesystem::path _get_exe() {

#ifdef _WIN32
    char path[MAX_PATH];
    GetModuleFileNameA(NULL, path, MAX_PATH);
    return std::filesystem::path(path);
#elif __APPLE__
    char path[1024];
    uint32_t size = sizeof(path);
    _NSGetExecutablePath(path, &size);
    return std::filesystem::path(path);
#elif __linux__
    char path[1024];
    ssize_t count = readlink("/proc/self/exe", path, sizeof(path));
    return std::filesystem::path(std::string(path, count));
#endif

}

#ifdef MPI_BUILD
#ifdef DEBUG
const std::string PROGRAM = "cmuts MPI (DEBUG)";
#else
const std::string PROGRAM = "cmuts MPI";
#endif
#else
#ifdef DEBUG
const std::string PROGRAM = "cmuts (DEBUG)";
#else
const std::string PROGRAM = "cmuts";
#endif
#endif

std::string _get_version() {

    auto exe = _get_exe().parent_path();
    auto version_file = exe / "VERSION";

    std::ifstream file(version_file);
    std::string version;
    std::getline(file, version);

    version.erase(version.find_last_not_of(" \t\n\r") + 1);
    version.erase(0, version.find_first_not_of(" \t\n\r"));

    return version;

}

void version() {
    std::cout << __repeat(SPACE, 8) + PROGRAM + " version " + _get_version() + "\n";
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
