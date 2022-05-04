#include <string>
#include <vector>

#include <h5.hpp>
#include <md.hpp>

#include "simulation_config.hpp"
#include "simulation_data.hpp"
#include "simulation_store.hpp"




simulation_store::simulation_store(std::string const& filename)
    : _file{filename, "w"}
{
}


void
simulation_store::save_config(simulation_config const& config)
{
    auto dataset = _file.dataset<h5::str>("metadata/config");
    dataset.write(dump_simulation_config(config));
}


void
simulation_store::save_beads(md::array_view<bead_data const> beads)
{
    auto dataset = _file.dataset<float, 2>("metadata/ab_factors");

    std::vector<float> values;
    for (auto const& bead : beads) {
        values.push_back(float(bead.a_factor));
        values.push_back(float(bead.b_factor));
    }
    dataset.write(values.data(), {beads.size(), 2});
}


void
simulation_store::save_chains(md::array_view<chain_data const> chains)
{
    auto dataset = _file.dataset<int, 2>("metadata/chain_ranges");

    std::vector<int> values;
    for (auto const& chain : chains) {
        values.push_back(int(chain.start));
        values.push_back(int(chain.end));
    }
    dataset.write(values.data(), {chains.size(), 2});
}


void
simulation_store::save_snapshot(md::step step, md::array_view<md::point const> positions)
{
    auto const key = std::to_string(step);

    auto positions_dataset = _file.dataset<float, 2>("snapshots/" + key + "/positions");
    positions_dataset.write(
        reinterpret_cast<md::scalar const*>(positions.data()),
        {positions.size(), 3},
        {.compression = 1, .scaleoffset = 3}
    );

    // Update keys
    auto keys_dataset = _file.dataset<h5::str, 1>("snapshots/.steps");
    std::vector<std::string> keys(keys_dataset.shape().size());
    if (keys.size() > 0) {
        keys_dataset.read(keys.data(), {keys.size()});
    }
    keys.push_back(key);
    keys_dataset.write(keys.data(), {keys.size()});
}
