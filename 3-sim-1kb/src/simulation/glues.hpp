#include <memory>

#include "config.hpp"
#include "glues/glue_simulator.hpp"


std::shared_ptr<glue_simulator> make_glue_simulator(simulation_config const& config);