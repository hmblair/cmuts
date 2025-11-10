#include "mutex.hpp"

namespace cmuts::mutex {

Mutex::~Mutex() {

    if (fd != -1) {
        close(fd);
        unlink(file.c_str());
    }

}


Mutex lock(const std::string& file) {

    Mutex mutex;
    mutex.file = file + ".lock";

    mutex.fd = open(mutex.file.c_str(), O_CREAT | O_EXCL | O_WRONLY, 0644);
    if (mutex.fd == -1) {
        throw std::runtime_error("Another process is already processing " + file);
    }

    if (flock(mutex.fd, LOCK_EX | LOCK_NB) == -1) {
        close(mutex.fd);
        unlink(mutex.file.c_str());
        mutex.fd = -1;
        throw std::runtime_error("Failed to acquire lock for " + file);
    }

    return mutex;

}


bool check(const std::string& file) {

    std::string lock = file + ".lock";

    int fd = open(lock.c_str(), O_RDONLY);
    if (fd == -1) return false;

    bool locked = (flock(fd, LOCK_EX | LOCK_NB) == -1);
    close(fd);

    return locked;

}


bool wait(const std::string& file) {

    bool locked = false;
    while (check(file)) {
        locked = true;
        usleep(100000);
    }
    return locked;

}

} // namespace cmuts::mutex
