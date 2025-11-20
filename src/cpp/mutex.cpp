#include "mutex.hpp"

namespace cmuts::mutex {

Mutex::~Mutex() {

    if (sem != SEM_FAILED) {
        sem_close(sem);
        sem_unlink(name.c_str());
    }

}

Mutex lock(const std::string& name) {

    Mutex mutex;
    mutex.name = "/" + name + "_lock";

    mutex.sem = sem_open(mutex.name.c_str(), O_CREAT | O_EXCL, 0644, 1);
    if (mutex.sem == SEM_FAILED) {
        throw std::runtime_error("Another process is already processing " + name);
    }

    if (sem_trywait(mutex.sem) == -1) {
        sem_close(mutex.sem);
        sem_unlink(mutex.name.c_str());
        mutex.sem = SEM_FAILED;
        throw std::runtime_error("Failed to acquire lock for " + name);
    }

    return mutex;

}

bool check(const std::string& name) {

    std::string semname = "/" + name + "_lock";

    sem_t* sem = sem_open(semname.c_str(), 0);
    if (sem == SEM_FAILED) return false;

    bool locked = (sem_trywait(sem) == -1);
    if (!locked) {
        sem_post(sem);
    }
    sem_close(sem);

    return locked;

}

bool wait(const std::string& name) {

    bool was_locked = false;
    while (check(name)) {
        was_locked = true;
        usleep(100000);
    }
    return was_locked;

}

} // namespace cmuts::mutex
