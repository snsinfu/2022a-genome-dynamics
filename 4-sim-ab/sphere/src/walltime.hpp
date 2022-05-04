#pragma once

#include <ctime>
#include <ostream>


// iostream-formattable wallclock time object.
class walltime
{
public:
    // Initializes object to the unix epoch.
    walltime() = default;

    // Initializes object to the given time.
    explicit
    walltime(std::time_t t);

    // Creates walltime object of current time,
    static
    walltime now();

    // Returns time_t value.
    std::time_t time() const;

private:
    std::time_t _time = 0;
};


std::ostream& operator<<(std::ostream& os, walltime const& t);
