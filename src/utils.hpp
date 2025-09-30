#ifndef CMUTS_UTIL_HEADER
#define CMUTS_UTIL_HEADER

#include <string>
#include <iostream>
#include <fstream>
#include <ctime>
#include <execinfo.h>
#include <argparse/argparse.hpp>

#ifdef _WIN32
    #include <windows.h>
#elif __APPLE__
    #include <mach-o/dyld.h>
#elif __linux__
    #include <unistd.h>
#endif

#include "mpi.hpp"

const std::string _LOG_FILE = "cmuts.log";
constexpr int32_t MAX_TRACE = 256;

typedef argparse::ArgumentParser Parser;

void __log(const std::string& filename, const std::string& message);
void __init_log(const std::string& filename);
void __throw_and_log(const std::string& filename, const std::string& err);

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
        T default_value
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
