#include <ctime>
#include <iomanip>
#include <ostream>

#include "walltime.hpp"


walltime::walltime(std::time_t t)
    : _time{t}
{
}

std::time_t
walltime::time() const
{
    return _time;
}


walltime
walltime::now()
{
    return walltime{std::time(nullptr)};
}


std::ostream&
operator<<(std::ostream& os, walltime const& t)
{
    auto unix_time = t.time();
    return os << std::put_time(std::localtime(&unix_time), "%F %T");
}
