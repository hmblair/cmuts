#include "mpi.hpp"

namespace MPI {

static inline int64_t ceil_div(int64_t x, int64_t y) {
    return x/y + (x % y != 0);
}
static inline int64_t round_up_to_multiple(int64_t x, int64_t y) {
    return ceil_div(x, y) * y;
}

std::string SPACE = " ";
std::string DIV   = "─";

static inline std::string __repeat(const std::string& str, size_t count) {

    std::string result;
    result.reserve(str.size() * count);
    for (size_t ix = 0; ix < count; ++ix) { result += str; }
    return result;

}

static inline std::string divider() {

    return __repeat(SPACE, 6) + __repeat(DIV, 35);

}


// Timer

static inline double __time() {
    #ifdef MPI_BUILD
    return MPI_Wtime();
    #else
    auto tp = std::chrono::high_resolution_clock::now();
    return std::chrono::duration<double>(tp.time_since_epoch()).count();
    #endif
}

static inline std::string _time_string(double elapsed) {

    int _seconds = static_cast<int>(elapsed);
    int hours = _seconds / HOUR;
    int minutes = (_seconds % HOUR) / MINUTE;
    int seconds = _seconds % MINUTE;

    std::stringstream ss;
    ss << std::setfill('0') << std::setw(2) << hours << ":"
       << std::setfill('0') << std::setw(2) << minutes << ":"
       << std::setfill('0') << std::setw(2) << seconds;
    return ss.str();

}

void Timer::start() {

    if (!running) {
        _start = __time();
        running = true;
    }

}

void Timer::stop() {

    if (running) {
        _total += __time() - _start;
        running = false;
    }

}

double Timer::elapsed() const {

    if (running) {
        return _total + (__time() - _start);
    }
    return _total;

}

std::string Timer::str() const {

    return _time_string(elapsed());

}

// Manager

#ifdef MPI_BUILD

static inline void __set_backing_dir() {
    setenv("OMPI_MCA_btl_sm_backing_directory", "/tmp", 0);
}

static inline void _init_thread(int& argc, char**& argv) {
    int _provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &_provided);
    if (_provided < MPI_THREAD_FUNNELED) {
        throw std::runtime_error("MPI error: Required threading level not supported.");
    }
}

#endif

Manager::Manager(int& argc, char**& argv) {

    #ifdef MPI_BUILD
    __set_backing_dir();
    _init_thread(argc, argv);
    MPI_Comm_rank(_comm, &_rank);
    MPI_Comm_size(_comm, &_size);
    MPI_Info_create(&_info);
    #endif

    timer.start();

}

Manager::~Manager() {

    timer.stop();

    #ifdef MPI_BUILD
    MPI_Finalize();
    #endif

}

int Manager::rank() const { return _rank; }

int Manager::size() const { return _size; }

#ifdef MPI_BUILD

MPI_Comm Manager::comm() const { return _comm; }

MPI_Info Manager::info() const { return _info; }

#endif

int64_t Manager::reduce(const int64_t& value) const {

    #ifdef MPI_BUILD

    int64_t _value = 0;
    if (!null()) {
        MPI_Reduce(&value, &_value, 1, MPI_INT, MPI_SUM, 0, _comm);
    }
    return _value;

    #else

    return value;

    #endif

}

bool Manager::root() const { return _rank == 0; }

bool Manager::null() const {

    #ifdef MPI_BUILD

    return _comm == MPI_COMM_NULL;

    #else

    return false;

    #endif

}

void Manager::barrier(bool debug) const {

    #ifdef MPI_BUILD

    if (debug) {
        std::cout << "Thread " << _rank << " waiting.\n";
    }

    if (!null()) { MPI_Barrier(_comm); }

    if (debug) {
        std::cout << "Thread " << _rank << " done.\n";
    }

    #else

    return;

    #endif

}

void Manager::remove(const std::string& name) const {

    if (root() && std::filesystem::is_regular_file(name)) {
        std::filesystem::remove(name);
    }
    barrier();

}

double Manager::time() const { return timer.elapsed(); }

std::string Manager::time_str() const { return timer.str(); }

void Manager::up(int lines) const {

    if (root()) {
        std::cout << "\033[" + std::to_string(lines) + "A";
    }

}

void Manager::down(int lines) const {

    if (root()) {
        for (int ix = 0; ix < lines; ix++) {
            std::cout << "\n";
        }
    }

}

void Manager::divide() const {

    if (root()) { std::cout << divider() << "\n"; }

}

int64_t Manager::chunksize(int64_t curr_size, int64_t total) const {

    return std::min(curr_size, total / _size);

}

Chunk Manager::chunk(int64_t _size, int64_t _total) const {

    Chunk chunk;

    chunk.low  = _size * rank();
    chunk.high = round_up_to_multiple(_total, _size * size());
    chunk.step = _size * size();
    chunk.size = _size;

    return chunk;

}

} // namespace MPI
