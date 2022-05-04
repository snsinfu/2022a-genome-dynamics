#pragma once

#include <memory>

#include "config.hpp"
#include "inits/initializer.hpp"


std::unique_ptr<simulation_initializer> make_initializer(simulation_config const& config);
