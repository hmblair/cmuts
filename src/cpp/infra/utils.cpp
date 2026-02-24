#include "infra/utils.hpp"

#include <sstream>

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

// Helper to call add_argument with or without a short name.
// When short_name is empty, passing it to argparse causes it to be
// stored in the names list, which breaks error messages (the empty
// string sorts first, producing "Error: : required." instead of
// "Error: --flag: required.").
static inline auto& _add(Parser& parser, const std::string& short_name, const std::string& long_name) {
    if (short_name.empty()) {
        return parser.add_argument(long_name);
    }
    return parser.add_argument(short_name, long_name);
}

template<>
void _add_arg<bool>(
    Parser& parser,
    const std::string& short_name,
    const std::string& long_name,
    const std::string& help
) {
    _add(parser, short_name, long_name)
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
    _add(parser, short_name, long_name)
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
    _add(parser, short_name, long_name)
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
    _add(parser, short_name, long_name)
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
    _add(parser, short_name, long_name)
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
    _add(parser, short_name, long_name)
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
    _add(parser, short_name, long_name)
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
    _add(parser, short_name, long_name)
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
    _add(parser, short_name, long_name)
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
    _add(parser, short_name, long_name)
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
