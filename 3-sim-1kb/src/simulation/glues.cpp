#include "glues.hpp"


std::shared_ptr<glue_simulator>
make_glue_simulator(simulation_config const& config)
{
    auto const box_size = config.chain.box_size;

    glue_simulator simulator({
        .max_glues      = config.glue.max_glues,
        .max_distance   = config.glue.glue_distance,
        .binding_rate   = config.glue.glue_binding_rate,
        .unbinding_rate = config.glue.glue_unbinding_rate,
        .box            = {.x_period = box_size, .y_period = box_size, .z_period = box_size},
    });

    return std::make_shared<glue_simulator>(simulator);
}