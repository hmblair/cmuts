#ifndef _CMUTS_MUTEX_HPP_
#define _CMUTS_MUTEX_HPP_

#include <string>
#include <stdexcept>
#include <fcntl.h>
#include <unistd.h>
#include <sys/file.h>

namespace cmuts::mutex {

struct Mutex {

    int fd = -1;
    std::string file;

    ~Mutex();

};


Mutex lock(const std::string& file);
bool check(const std::string& file);
bool wait(const std::string& file);

} // namespace cmuts::mutex

#endif
