#pragma once

#include <cstddef>
#include <random>
#include <vector>

#include "loop_simulator.hpp"


class basic_loop_simulator : public loop_simulator
{
public:
    struct constructor_config
    {
        /** The length of a polymer chain which loop factors move along. */
        std::size_t chain_length = 0;

        /** The maximum number of simulated loops. */
        std::size_t max_loops = 0;
    };

    explicit basic_loop_simulator(constructor_config const& config);

    /** Sets kinetic loading rate of loop factors. */
    void set_loading_rate(double r);

    /** Sets kinetic unloading rate of loop factors. */
    void set_unloading_rate(double r);

    /** Sets speed of loop factors in the extending direction. */
    void set_forward_speed(double v);

    /** Sets speed of loop factors in the shrinking direction. */
    void set_backward_speed(double v);

    /** Sets kinetic crossing (Z-looping) rate of loop factors upon collision. */
    void set_crossing_rate(double r);

    /** Sets relative speed of loop factors entering to specified position. */
    void set_site_attachability(std::size_t pos, double m);

    /** Sets relative speed of loop factors leaving from specified position. */
    void set_site_detachability(std::size_t pos, double m);

    /** Adds boundary at the specified position. */
    void add_boundary(std::size_t pos);

    /** Forcifully loads a new handcuff loop at specified position. */
    void load_loop(std::size_t pos);

    std::size_t      chain_length() const;
    void             clear() override;
    loop_pair const* begin() const override;
    loop_pair const* end() const override;
    void             step(double dt, std::mt19937_64& random) override;
    void             preload(std::mt19937_64& random) override;

private:
    void step_unloading(double dt, std::mt19937_64& random);
    void step_loading(double dt, std::mt19937_64& random);
    void step_single_loading(double dt, std::mt19937_64& random);
    void step_motion(double dt, std::mt19937_64& random);

private:
    std::vector<loop_pair> _loops;
    std::vector<double>    _sites_attachability;
    std::vector<double>    _sites_detachability;
    std::vector<unsigned>  _sites_state;
    std::vector<int>       _sites_occupancy;
    std::size_t            _next_id = 1;
    double                 _loading_rate = 0;
    double                 _unloading_rate = 0;
    double                 _forward_speed = 0;
    double                 _backward_speed = 0;
    double                 _crossing_rate = 0;
    bool                   _check_crossing = false;
};
