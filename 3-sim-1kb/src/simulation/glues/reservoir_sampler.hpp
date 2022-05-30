#pragma once

#include <cmath>
#include <cstddef>
#include <random>
#include <vector>


/**
 * Container-like class for collecting fixed number of samples uniformly from a
 * sequence of unknown length.
 */
template<typename T>
class reservoir_sampler
{
public:
    using iterator = T const*;

    /** Creates a sampler with given capacity. */
    explicit reservoir_sampler(std::size_t capacity)
        : _capacity{capacity}
    {
        _samples.reserve(capacity);
    }

    /** Returns the maximum number of samples the sampler holds. */
    std::size_t capacity() const
    {
        return _capacity;
    }

    /** Returns the number of samples the sampler currently holds. */
    std::size_t size() const
    {
        return _samples.size();
    }

    /** Returns the number of samples the sampler has seen. */
    std::size_t population_size() const
    {
        return _population_size;
    }

    /** Returns the beginning of the sequence of samples held by the sampler. */
    iterator begin() const
    {
        return _samples.data();
    }

    /** Returns the past-the-end of the sequence of samples held by the sampler. */
    iterator end() const
    {
        return _samples.data() + _samples.size();
    }

    /** Feeds a single sample. */
    template<typename RNG>
    void feed(T const& value, RNG& random)
    {
        // Implements reservoir sampling "Algorithm-L".

        _population_size++;

        if (_population_size <= _capacity) [[unlikely]] {
            _samples.push_back(value);
            return;
        }

        if (_population_size == _capacity + 1) [[unlikely]] {
            update_skip(random);
        }

        if (_skip != 0) [[likely]] {
            _skip--;
            return;
        }

        std::uniform_int_distribution<std::size_t> slot{0, _capacity - 1};
        _samples[slot(random)] = value;

        update_skip(random);
    }

private:
    template<typename RNG>
    void update_skip(RNG& random)
    {
         std::uniform_real_distribution<double> uniform;
         _acceptance *= std::pow(uniform(random), 1 / double(_capacity));

         std::geometric_distribution<std::size_t> skip{_acceptance};
         _skip = skip(random);
    }

private:
    std::size_t    _capacity;
    std::vector<T> _samples;
    std::size_t    _population_size = 0;
    std::size_t    _skip = 0;
    double         _acceptance = 1;
};
