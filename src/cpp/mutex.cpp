#include "mutex.hpp"

namespace cmuts::mutex {

Mutex::~Mutex() {
    if (fd != -1) {
        close(fd);
    }
}

Mutex lock(const std::string& name) {

    Mutex mutex;
    mutex.name = name;

    mutex.fd = open(name.c_str(), O_RDWR);
    if (mutex.fd == -1) {
        throw std::runtime_error("Failed to open " + name + " for locking");
    }

    if (flock(mutex.fd, LOCK_EX | LOCK_NB) == -1) {
        close(mutex.fd);
        mutex.fd = -1;
        throw std::runtime_error("Another process is already processing " + name);
    }

    return mutex;

}

bool check(const std::string& name) {

    int fd = open(name.c_str(), O_RDONLY);
    if (fd == -1) return false;

    bool locked = (flock(fd, LOCK_EX | LOCK_NB) == -1);
    if (!locked) {
        flock(fd, LOCK_UN);
    }
    close(fd);

    return locked;

}

bool wait(const std::string& name) {

    int fd = open(name.c_str(), O_RDONLY);
    if (fd == -1) return false;

    // Block until we can acquire the lock, then release immediately
    bool was_locked = (flock(fd, LOCK_EX) == 0);
    if (was_locked) {
        flock(fd, LOCK_UN);
    }
    close(fd);

    return was_locked;

}

} // namespace cmuts::mutex
