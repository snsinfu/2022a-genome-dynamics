#pragma once

#include <memory>

#include "config.hpp"
#include "loops/loop_simulator.hpp"


std::shared_ptr<loop_simulator> make_loop_simulator(simulation_config const& config);
