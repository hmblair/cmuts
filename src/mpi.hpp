#ifndef MPI_HEADER
#define MPI_HEADER

#include <filesystem>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#ifdef MPI_BUILD
extern "C" {
    #include <mpi.h>
}
#else
#include <chrono>
#endif

const int MINUTE = 60;
const int HOUR   = 3600;

namespace MPI {

class Timer {
private:

    bool running  = false;
    double _start = 0.0;
    double _total = 0.0;

public:

    Timer() = default;
    void start();
    void stop();
    double elapsed() const;
    std::string str() const;

};

struct Chunk {

    int64_t low  = 0;
    int64_t high = 0;
    int64_t step = 0;
    int64_t size = 0;

};

class Manager {
private:

    int _rank = 0;
    int _size = 1;

    #ifdef MPI_BUILD
    MPI_Info _info = nullptr;
    MPI_Comm _comm = MPI_COMM_WORLD;
    #endif

    Timer timer;

public:

    Manager(int& argc, char**& argv);
    Manager (const Manager&) = delete;
    Manager& operator= (const Manager&) = delete;
    ~Manager();

    int rank() const;
    int size() const;

    #ifdef MPI_BUILD
    MPI_Comm comm() const;
    MPI_Info info() const;
    #endif

    int64_t reduce(const int64_t& value) const;
    int64_t broadcast(int64_t& value) const;
    bool root() const;
    bool null() const;
    void barrier(bool debug = false) const;
    void remove(const std::string& file) const;
    void print(const std::string& str) const;
    double time() const;
    std::string time_str() const;
    void up(int lines = 1) const;
    void down(int lines = 1) const;
    void divide() const;

    class OutStream {
    public:
        OutStream(const Manager& manager) : manager(manager) {}
        ~OutStream() {
            if (manager.root() && !oss.str().empty()) {
                std::cout << oss.str();
            }
        }

        template <typename T>
        OutStream& operator<<(const T& val) {
            oss << val;
            return *this;
        }

    private:
        const Manager& manager;
        std::ostringstream oss;
    };
    OutStream out() const {
        return OutStream(*this);
    }

    class ErrStream {
    public:
        ErrStream(const Manager& manager) : manager(manager) {}
        ~ErrStream() {
            if (manager.root() && !oss.str().empty()) {
                std::cerr << oss.str();
            }
        }

        template <typename T>
        ErrStream& operator<<(const T& val) {
            oss << val;
            return *this;
        }

    private:
        const Manager& manager;
        std::ostringstream oss;
    };
    ErrStream err() const {
        return ErrStream(*this);
    }

    Chunk chunk(int64_t size, int64_t total) const;

};

} // namespace MPI

#endif
