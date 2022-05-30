#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <set>
#include <string>
#include <vector>

#include <highfive/H5Attribute.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>
#include <highfive/H5PropertyList.hpp>
#include <nlohmann/json.hpp>
#include <md.hpp>

#include "h5/vector_array.hpp"
#include "particle_data.hpp"
#include "simulation_store.hpp"


namespace
{
    template<typename T, std::size_t cols>
    H5::DataSet write_compressed_array(
        H5::Group& group,
        std::string const& name,
        std::vector<std::array<T, cols>> const& data
    );

    void clear_dataset(H5::Group& group, std::string const& name);

    template<typename G>
    H5::Group require_group(G& parent, std::string name);

    template<typename G>
    H5::Group require_snapshot_group(G& root, std::string phase, md::step step);

    float quantize(md::scalar val, int bits);
}


simulation_store::simulation_store(std::string const& filename)
    : _store{filename, H5::File::ReadWrite}
{
}


simulation_config simulation_store::load_config()
{
    auto metadata = _store.getGroup("metadata");
    auto config_data = metadata.getDataSet("config");

    std::string config_str;
    config_data.read(config_str);

    simulation_config config;
    parse_simulation_config(config_str, config);
    return config;
}


std::vector<chromosome_range> simulation_store::load_chromosomes()
{
    auto metadata = _store.getGroup("metadata");
    auto chromosome_ranges_data = metadata.getDataSet("chromosome_ranges");
    auto centromere_ranges_data = metadata.getDataSet("centromere_ranges");

    std::string keys_json;
    chromosome_ranges_data.getAttribute("keys").read(keys_json);
    auto const keys = nlohmann::json::parse(keys_json);

    std::vector<std::array<int, 2>> chromosome_ranges;
    std::vector<std::array<int, 2>> centromere_ranges;
    chromosome_ranges_data.read(chromosome_ranges);
    centromere_ranges_data.read(centromere_ranges);

    std::vector<chromosome_range> chroms;
    chroms.reserve(chromosome_ranges.size());

    for (std::size_t i = 0; i < chromosome_ranges.size(); i++) {
        std::string name;
        for (auto it = keys.begin(); it != keys.end(); it++) {
            if (it.value() == i) {
                name = it.key();
                break;
            }
        }

        chroms.push_back({
            .name             = name,
            .start            = static_cast<md::index>(chromosome_ranges[i][0]),
            .end              = static_cast<md::index>(chromosome_ranges[i][1]),
            .centromere_start = static_cast<md::index>(centromere_ranges[i][0]),
            .centromere_end   = static_cast<md::index>(centromere_ranges[i][1])
        });
    }

    return chroms;
}


std::vector<particle_data> simulation_store::load_particle_data()
{
    auto metadata = _store.getGroup("metadata");
    auto ab_factors_data = metadata.getDataSet("ab_factors");

    std::vector<std::array<md::scalar, 2>> ab_factors;
    ab_factors_data.read(ab_factors);

    std::vector<particle_data> particles;
    particles.reserve(ab_factors.size());

    for (auto const& ab : ab_factors) {
        particles.push_back({
            .a_factor = ab[0],
            .b_factor = ab[1]
        });
    }

    return particles;
}


std::vector<index_range> simulation_store::load_nucleolus_ranges()
{
    std::vector<std::array<int, 2>> range_values;
    auto metadata = _store.getGroup("metadata");
    auto nucleolus_ranges_data = metadata.getDataSet("nucleolus_ranges");
    if (nucleolus_ranges_data.getElementCount() > 0) {
        // Some simulation assumes a system with no nucleolar particle. In
        // that case the dataspace of nucleolus_ranges dataset is {0, 0}.
        // HighFive throws an exception on reading such a dataset, so handle
        // empty dataset ourselves.
        nucleolus_ranges_data.read(range_values);
    }

    std::vector<index_range> ranges;
    for (auto const& [start, end] : range_values) {
        ranges.push_back({
            .begin = static_cast<md::index>(start),
            .end   = static_cast<md::index>(end)
        });
    }

    return ranges;
}


std::vector<nucleolus_bond> simulation_store::load_nucleolus_bonds()
{
    std::vector<std::array<int, 2>> index_pairs;
    auto metadata = _store.getGroup("metadata");
    auto nucleolus_bond_data = metadata.getDataSet("nucleolus_bonds");
    if (nucleolus_bond_data.getElementCount() > 0) {
        nucleolus_bond_data.read(index_pairs);
    }

    std::vector<nucleolus_bond> bonds;
    for (auto const& pair : index_pairs) {
        bonds.push_back({
            .nor_index = static_cast<md::index>(pair[0]),
            .nuc_index = static_cast<md::index>(pair[1])
        });
    }

    return bonds;
}


void simulation_store::set_phase(std::string const& phase)
{
    _phase = phase;
}


void simulation_store::save_chromosomes(md::array_view<chromosome_range const> chroms)
{
    auto snapshots_group = require_group(_store, "snapshots");
    auto phase_group = require_group(snapshots_group, _phase);
    auto metadata_group = require_group(phase_group, "metadata");

    std::vector<std::array<int, 2>> ranges;
    nlohmann::json keys;

    for (auto const& chrom : chroms) {
        ranges.push_back({
            static_cast<int>(chrom.start),
            static_cast<int>(chrom.end)
        });
        keys[chrom.name] = ranges.size() - 1; // index
    }

    clear_dataset(metadata_group, "chromosome_ranges");
    auto dataset = write_compressed_array(metadata_group, "chromosome_ranges", ranges);
    dataset.createAttribute("keys", keys.dump());

    _store.flush();
}


void simulation_store::save_context(md::step step, simulation_context const& context)
{
    auto snapshot = require_snapshot_group(_store, _phase, step);

    std::vector<md::scalar> const wall_semiaxes = {
        context.wall_semiaxes.x,
        context.wall_semiaxes.y,
        context.wall_semiaxes.z
    };

    nlohmann::json json;
    json["time"]          = context.time;
    json["bead_scale"]    = context.bead_scale;
    json["bond_scale"]    = context.bond_scale;
    json["wall_semiaxes"] = wall_semiaxes;
    json["mean_energy"]   = context.mean_energy;
    json["wall_energy"]   = context.wall_energy;

    clear_dataset(snapshot, "context");
    snapshot.createDataSet("context", json.dump());

    // Commit changes to the disk so that the HDF5 file can be inspected before
    // the simulation finishes.
    _store.flush();
}


simulation_context simulation_store::load_context(md::step step)
{
    auto snapshot = require_snapshot_group(_store, _phase, step);

    std::string context_json;
    snapshot.getDataSet("context").read(context_json);

    auto json = nlohmann::json::parse(context_json.begin(), context_json.end());

    simulation_context context;
    context.time        = json["time"];
    context.bead_scale  = json["bead_scale"];
    context.bond_scale  = json["bond_scale"];
    context.mean_energy = json["mean_energy"];
    context.wall_energy = json["wall_energy"];

    std::vector<md::scalar> vec = json["wall_semiaxes"];
    context.wall_semiaxes = {vec[0], vec[1], vec[2]};

    return context;
}


void simulation_store::save_positions(md::step step, md::array_view<md::point const> positions)
{
    auto snapshot = require_snapshot_group(_store, _phase, step);

    // Quantize coordinate values for better compression. Resolution of
    // ~0.0001 (0.1 nm) is sufficient for our simulation, so use 16 bits.
    constexpr int fraction_bits = 16;

    std::vector<std::array<float, 3>> positions_array(positions.size());
    for (md::index i = 0; i < positions.size(); i++) {
        positions_array[i] = {
            quantize(positions[i].x, fraction_bits),
            quantize(positions[i].y, fraction_bits),
            quantize(positions[i].z, fraction_bits)
        };
    }

    clear_dataset(snapshot, "positions");
    write_compressed_array(snapshot, "positions", positions_array);

    // Commit changes to the disk so that the HDF5 file can be inspected before
    // the simulation finishes.
    _store.flush();
}


std::vector<md::point> simulation_store::load_positions(md::step step)
{
    auto snapshot = require_snapshot_group(_store, _phase, step);

    std::vector<std::array<float, 3>> positions_array;
    snapshot.getDataSet("positions").read(positions_array);

    std::vector<md::point> positions;
    positions.reserve(positions_array.size());
    for (auto const& coords : positions_array) {
        positions.push_back({
            static_cast<md::scalar>(coords[0]),
            static_cast<md::scalar>(coords[1]),
            static_cast<md::scalar>(coords[2])
        });
    }

    return positions;
}


void simulation_store::save_contacts(
    md::step step, std::vector<std::array<std::uint32_t, 3>> const& contacts
)
{
    if (contacts.empty()) {
        return;
    }

    auto snapshot = require_snapshot_group(_store, _phase, step);
    clear_dataset(snapshot, "contact_map");
    write_compressed_array(snapshot, "contact_map", contacts);

    _store.flush();
}


namespace
{
    // Aim for 1MB.
    constexpr std::size_t h5_chunk_size = 1024 * 1024;


    template<typename T, std::size_t cols>
    H5::DataSet write_compressed_array(
        H5::Group& group,
        std::string const& name,
        std::vector<std::array<T, cols>> const& data
    )
    {
        auto const chunk_rows = std::min(
            h5_chunk_size / sizeof(T[cols]),
            data.size()
        );

        H5::DataSpace dataspace(data.size(), cols);
        H5::DataSetCreateProps props;
        props.add(H5::Chunking(chunk_rows, cols));
        props.add(H5::Shuffle());
        props.add(H5::Deflate(6));

        auto dataset = group.createDataSet<T>(name, dataspace, props);
        dataset.write(data);

        return dataset;
    }


    void clear_dataset(H5::Group& group, std::string const& name)
    {
        if (group.exist(name)) {
            (void) H5Gunlink(group.getId(), name.c_str());
        }
    }


    template<typename G>
    H5::Group require_group(G& parent, std::string name)
    {
        if (!parent.exist(name)) {
            return parent.createGroup(name);
        }
        return parent.getGroup(name);
    }


    template<typename G>
    void update_ordered_steps(G& group, md::step new_step)
    {
        constexpr char const* dataset_name = ".steps";

        std::vector<std::string> steps;
        if (group.exist(dataset_name)) {
            group.getDataSet(dataset_name).read(steps);
        } else {
            group.createDataSet(dataset_name, steps);
        }

        std::set<long> sorted_steps;
        for (auto const step : steps) {
            sorted_steps.insert(std::stol(step));
        }
        sorted_steps.insert(static_cast<long>(new_step));

        steps.clear();
        for (auto const step : sorted_steps) {
            steps.push_back(std::to_string(step));
        }

        clear_dataset(group, dataset_name);
        group.createDataSet(dataset_name, steps);
    }

    template<typename G>
    H5::Group require_snapshot_group(G& root, std::string phase, md::step step)
    {
        auto snapshots_group = require_group(root, "snapshots");
        auto phase_group = require_group(snapshots_group, phase);
        auto snapshot_group = require_group(phase_group, std::to_string(step));
        update_ordered_steps(phase_group, step);
        return snapshot_group;
    }


    float quantize(md::scalar val, int bits)
    {
        float const scale = 1 << bits;
        return std::nearbyint(float(val) * scale) / scale;
    }
}
