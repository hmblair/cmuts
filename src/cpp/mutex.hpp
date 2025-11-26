#ifndef _CMUTS_MUTEX_HPP_
#define _CMUTS_MUTEX_HPP_

#include <string>
#include <stdexcept>
#include <unistd.h>
#include <errno.h>
#include <fcntl.h>
#include <sys/file.h>

namespace cmuts::mutex {

struct Mutex {
    int fd = -1;
    std::string name;
    ~Mutex();
};

Mutex lock(const std::string& file);
bool check(const std::string& file);
bool wait(const std::string& file);

} // namespace cmuts::mutex

#endif
