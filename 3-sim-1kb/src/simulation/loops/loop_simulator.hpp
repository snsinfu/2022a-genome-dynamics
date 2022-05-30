#pragma once

#include <cstddef>
#include <random>


/**
 * State, or the pair of bound positions and the ID number, of a single cohesin
 * loop. Each active loop has unique ID number >= 1.
 */
struct loop_pair
{
    std::size_t start = 0;
    std::size_t end   = 0;
    std::size_t id    = 0;
};


/** Generic interface of loop simulator. */
class loop_simulator
{
public:
    virtual ~loop_simulator() = default;

    /** Unloads all loops. */
    virtual void clear() = 0;

    /** Returns the start of the sequence of simulated loops. */
    virtual loop_pair const* begin() const = 0;

    /** Returns past the end of the sequence of simulated loops. */
    virtual loop_pair const* end() const = 0;

    /** Returns the number of simulated loops. */
    std::size_t size() const
    {
        return static_cast<std::size_t>(end() - begin());
    }

    /** Simulates the motion of the loops within given time interval. */
    virtual void step(double dt, std::mt19937_64& random) = 0;

    /** Attempts to load specified number of loops at random positions. */
    virtual void preload(std::mt19937_64& random) = 0;
};
