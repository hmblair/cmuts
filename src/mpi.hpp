#ifndef MPI_HEADER
#define MPI_HEADER

#include <stdexcept>
#include <string>
#include <sstream>
#include <iostream>
#include <filesystem>
extern "C" {
    #include <mpi.h>
}

const int MS     = 1000;
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
    int _size = 0;
    MPI_Info _info = nullptr;
    MPI_Comm _comm = MPI_COMM_WORLD;
    Timer timer;

public:

    Manager(int& argc, char**& argv);
    Manager (const Manager&) = delete;
    Manager& operator= (const Manager&) = delete;
    ~Manager();

    int rank() const;
    int size() const;
    MPI_Comm comm() const;
    MPI_Info info() const;

    int64_t reduce(const int64_t& value) const;
    int64_t broadcast(int64_t& value) const;
    bool root() const;
    bool null() const;
    void barrier() const;
    void remove(const std::string& file) const;
    void print(const std::string& str) const;
    double time() const;
    std::string time_str() const;

    class ErrStream {
    public:
        ErrStream(const Manager& manager) : manager(manager) {}
        ~ErrStream() {
            if (manager.root() && !oss.str().empty()) {
                std::cerr << oss.str();
            }
            manager.barrier();
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
