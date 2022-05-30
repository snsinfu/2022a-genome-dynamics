#include <string>

#include <h5.hpp>

#include "buffer_traits.hpp"
#include "config.hpp"
#include "store.hpp"
#include "topology.hpp"


simulation_store::simulation_store(std::string const& filename)
    : _file{filename, "w"}
{
}


void
simulation_store::save_metadata(metadata_record const& metadata)
{
    _file.dataset<h5::str>("config").write(
        format_simulation_config(metadata.config)
    );
    _file.dataset<h5::str>("config_source").write(metadata.config.config_text);

    std::vector<md::index> flat_chain_ranges;
    for (auto const& chain : make_chain_assignments(metadata.config)) {
        flat_chain_ranges.push_back(chain.start);
        flat_chain_ranges.push_back(chain.end);
    }
    _file.dataset<h5::i32, 2>("chain_ranges").write(
        flat_chain_ranges.data(), {flat_chain_ranges.size() / 2, 2}
    );
}


void
simulation_store::save_snapshot(snapshot_record const& snapshot)
{
    if (!_positions_dataset) {
        _positions_dataset = _file.dataset<h5::f32, 3>("positions_history");
        _positions_stream = _positions_dataset->stream_writer(
            h5::shape<2>{snapshot.positions.size(), 3}, {.compression = 1}
        );
    }

    if (!_loops_dataset && snapshot.loops.data()) {
        _loops_dataset = _file.dataset<h5::i32, 3>("loops_history");
        _loops_stream = _loops_dataset->stream_writer(
            h5::shape<2>{snapshot.loops.size(), 3}, {.compression = 1}
        );
    }

    _positions_stream->write(snapshot.positions);

    if (_loops_stream) {
        _loops_stream->write(snapshot.loops);
    }
}
