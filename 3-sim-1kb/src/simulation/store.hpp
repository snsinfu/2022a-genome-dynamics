#pragma once

#include <optional>
#include <string>
#include <vector>

#include <h5.hpp>
#include <md.hpp>

#include "config.hpp"
#include "loops/loop_simulator.hpp"


struct metadata_record
{
    simulation_config const& config;
};


struct snapshot_record
{
    md::array_view<md::point const> positions;
    md::array_view<loop_pair const> loops;
};


class simulation_store
{
public:
    explicit simulation_store(std::string const& filename);
    void     save_metadata(metadata_record const& metadata);
    void     save_snapshot(snapshot_record const& snapshot);

private:
    h5::file                                     _file;
    std::optional<h5::dataset      <h5::f32, 3>> _positions_dataset;
    std::optional<h5::stream_writer<h5::f32, 2>> _positions_stream;
    std::optional<h5::dataset      <h5::i32, 3>> _loops_dataset;
    std::optional<h5::stream_writer<h5::i32, 2>> _loops_stream;
};
