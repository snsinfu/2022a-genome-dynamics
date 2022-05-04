#include <algorithm>
#include <cmath>
#include <cstddef>
#include <random>

#include "basic_loop_simulator.hpp"


namespace site_state
{
    enum : unsigned
    {
        none     = 0,
        blocked  = 1 << 0,  // loop factors can never step into this site
        boundary = 1 << 1,  // site bears a boundary element
    };
}


basic_loop_simulator::basic_loop_simulator(constructor_config const& config)
    : _loops(config.max_loops)
    , _sites_attachability(config.chain_length, 1)
    , _sites_detachability(config.chain_length, 1)
    , _sites_state(config.chain_length, 0)
    , _sites_occupancy(config.chain_length, 0)
{
    // We have `max_loops` loops that are not yet loaded. Set dummy states.
    loop_pair const dummy = {
        .start = config.chain_length,
        .end   = config.chain_length,
        .id    = 0,
    };
    std::fill(_loops.begin(), _loops.end(), dummy);
}


void
basic_loop_simulator::set_loading_rate(double r)
{
    _loading_rate = r;
}


void
basic_loop_simulator::set_unloading_rate(double r)
{
    _unloading_rate = r;
}


void
basic_loop_simulator::set_forward_speed(double v)
{
    _forward_speed = v;
}


void
basic_loop_simulator::set_backward_speed(double v)
{
    _backward_speed = v;
}


void
basic_loop_simulator::set_crossing_rate(double r)
{
    // Inifinite crossing rate mathematically means that crossing must always
    // succeed. We special-case that as freely diffusing loop factors.
    if (std::isinf(r)) {
        _crossing_rate = 0;
        _check_crossing = false;
    } else {
        _crossing_rate = r;
        _check_crossing = true;
    }
}


void
basic_loop_simulator::set_site_attachability(std::size_t pos, double m)
{
    _sites_attachability[pos] = m;
}


void
basic_loop_simulator::set_site_detachability(std::size_t pos, double m)
{
    _sites_detachability[pos] = m;
}


void
basic_loop_simulator::add_boundary(std::size_t pos)
{
    _sites_state[pos] |= site_state::blocked | site_state::boundary;
}


void
basic_loop_simulator::load_loop(std::size_t pos)
{
    // Do not load a loop on boundary elements.
    if (_sites_state[pos] & site_state::blocked) {
        return;
    }

    // Find a usable loop instance to load. Linear search is ok because loop
    // loading is relatively infrequent.
    auto const unused_loop = std::find_if(
        _loops.begin(), _loops.end(), [](auto& loop) { return loop.id == 0; }
    );
    if (unused_loop == _loops.end()) {
        return;
    }

    // Load a loop factor. This function does not check site occupancy. The
    // caller is responsible for determining if loop overloading is possible.
    auto& loop = *unused_loop;
    loop = {
        .start = pos,
        .end   = pos,
        .id    = _next_id++,
    };
    _sites_occupancy[pos] += 2;
}


std::size_t
basic_loop_simulator::chain_length() const
{
    return _sites_state.size();
}


void
basic_loop_simulator::clear()
{
    loop_pair const dummy = {
        .start = chain_length(),
        .end   = chain_length(),
        .id    = 0,
    };

    std::fill(_loops.begin(), _loops.end(), dummy);
    std::fill(_sites_occupancy.begin(), _sites_occupancy.end(), 0);
}


loop_pair const*
basic_loop_simulator::begin() const
{
    return _loops.data();
}


loop_pair const*
basic_loop_simulator::end() const
{
    return _loops.data() + _loops.size();
}


void
basic_loop_simulator::step(double dt, std::mt19937_64& random)
{
    step_unloading(dt, random);
    step_loading(dt, random);
    step_motion(dt, random);
}


/** Simulates random unloading of loaded loop factors. */
void
basic_loop_simulator::step_unloading(double dt, std::mt19937_64& random)
{
    if (_unloading_rate == 0) {
        return;
    }

    loop_pair const dummy = {
        .start = chain_length(),
        .end   = chain_length(),
        .id    = 0,
    };

    for (auto& loop : _loops) {
        if (!loop.id) {
            continue;
        }

        // Unloading should slow down if cohesin is physically adsorbed on the
        // currently bound site.
        auto const rate = _unloading_rate * std::min(
            _sites_detachability[loop.start],
            _sites_detachability[loop.end]
        );
        std::bernoulli_distribution unloading_event{-std::expm1(-rate * dt)};

        if (unloading_event(random)) {
            _sites_occupancy[loop.start] -= 1;
            _sites_occupancy[loop.end] -= 1;
            loop = dummy;
        }
    }
}


/** Simulates random loading of new loop factors. */
void
basic_loop_simulator::step_loading(double dt, std::mt19937_64& random)
{
    if (_loading_rate == 0) {
        return;
    }

    std::poisson_distribution<int> loading_events{_loading_rate * dt};
    auto const count = loading_events(random);

    std::uniform_int_distribution<std::size_t> loading_site{0, chain_length() - 1};
    std::poisson_distribution<int> crossing{_crossing_rate * dt};

    for (int i = 0; i < count; i++) {
        auto const pos = loading_site(random);

        if (_sites_state[pos] & site_state::blocked) {
            continue;
        }

        // Treat overloading as collisions (loop crossings). So, overloading
        // successds at a predefined rate.
        if (_check_crossing && crossing(random) < _sites_occupancy[pos]) {
            continue;
        }

        std::bernoulli_distribution attaching{_sites_attachability[pos]};
        if (!attaching(random)) {
            continue;
        }

        load_loop(pos);
    }
}


/** Simulates random motion of loaded loop factors. */
void
basic_loop_simulator::step_motion(double dt, std::mt19937_64& random)
{
    auto const move_onto = [&](std::size_t& pos, std::size_t dest) {
        _sites_occupancy[pos] -= 1;
        _sites_occupancy[dest] += 1;
        pos = dest;
    };

    auto const attempt_move = [&](std::size_t& pos, std::size_t dest, double rate) {
        if (_sites_state[dest] & site_state::blocked) {
            return;
        }

        // Detachability affects all motion kinetics, so let us scale dt.
        auto const scaled_dt = dt * _sites_detachability[pos] * _sites_attachability[dest];

        // Handle two physically exclusive events: Crossing and free slip.
        std::poisson_distribution<int> crossing{_crossing_rate * scaled_dt};
        std::bernoulli_distribution slipping{-std::expm1(-rate * scaled_dt)};

        if (_check_crossing && _sites_occupancy[dest] > 0) {
            // Crossing occurs when loop roots collide. There may be multiple
            // (say n) roots occupying a site, so for another root to step onto
            // the site, crossing event must occur n times in a timestep.
            if (crossing(random) < _sites_occupancy[dest]) {
                return;
            }
        } else {
            // Kinetic slip onto an unoccupied site.
            if (!slipping(random)) {
                return;
            }
        }

        move_onto(pos, dest);
    };

    auto const attempt_move_left = [&](std::size_t& pos, double rate) {
        if (pos > 0) {
            attempt_move(pos, pos - 1, rate);
        }
    };

    auto const attempt_move_right = [&](std::size_t& pos, double rate) {
        if (pos + 1 < chain_length()) {
            attempt_move(pos, pos + 1, rate);
        }
    };

    for (auto& loop : _loops) {
        if (!loop.id) {
            continue;
        }

        // Kinetically move start (upstream) and end (downstream) roots in
        // their "forward" and "backward" directions.
        attempt_move_left(loop.start, _forward_speed);
        attempt_move_right(loop.start, _backward_speed);
        attempt_move_right(loop.end, _forward_speed);
        attempt_move_left(loop.end, _backward_speed);

        // Start and end roots bounce off each other in case of collision.
        if (loop.start > loop.end) {
            std::swap(loop.start, loop.end);
        }
    }
}


void
basic_loop_simulator::preload(std::mt19937_64& random)
{
    // Load expected number of loops randomly. Overloading is allowed.
    auto const expected_count = std::size_t(_loading_rate / _unloading_rate);

    std::uniform_int_distribution<std::size_t> loading_site{0, chain_length() - 1};

    for (std::size_t i = 0; i < expected_count; i++) {
        auto const pos = loading_site(random);

        if (_sites_state[pos] & site_state::blocked) {
            continue;
        }

        std::bernoulli_distribution attaching{_sites_attachability[pos]};
        if (!attaching(random)) {
            continue;
        }

        load_loop(pos);
    }
}
