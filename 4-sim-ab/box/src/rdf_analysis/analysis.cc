#include <cmath>
#include <cstddef>
#include <iostream>
#include <string>
#include <vector>

#include <h5.hpp>
#include <md.hpp>
#include <nlohmann/json.hpp>

#include "analysis.hpp"
#include "distance_histogram.hpp"


void run_analysis(analysis_config const& config)
{
    h5::file store{config.filename, "r"};

    // Load particle A/B parameters. The parameters are used to select
    // particles to analyze.
    std::vector<std::pair<md::scalar, md::scalar>> ab_factors;
    {
        auto dataset = store.dataset<float, 2>("metadata/ab_factors");
        ab_factors.resize(dataset.shape().dims[0]);
        dataset.read(
            reinterpret_cast<md::scalar*>(ab_factors.data()),
            {ab_factors.size(), 2}
        );
    }

    // Select particles. We store an array of the indices of selected particles.
    md::scalar a_factor = -1; // -1 means all

    if (config.particle_selection == "A") {
        a_factor = 1;
    }
    if (config.particle_selection == "B") {
        a_factor = 0;
    }

    std::vector<md::index> particle_selection;
    if (a_factor == -1) {
        for (md::index i = 0; i < ab_factors.size(); i++) {
            particle_selection.push_back(i);
        }
    } else {
        for (md::index i = 0; i < ab_factors.size(); i++) {
            if (std::fabs(ab_factors[i].first - a_factor) < 0.1) {
                particle_selection.push_back(i);
            }
        }
    }

    auto const num_points = particle_selection.size();

    // Calculate expected density of selected points. This is used as the prior
    // density for RDF computation.
    md::scalar box_size;
    {
        std::string config_json;
        store.dataset<h5::str>("metadata/config").read(config_json);
        auto const simulation_config = nlohmann::json::parse(config_json);
        box_size = simulation_config["box_size"];
    }
    auto const volume = md::power<3>(box_size);
    auto const expected_density = md::scalar(num_points) / volume;

    // Select snapshots to analyze. TODO: Use config.step_start/step_end.
    std::vector<std::string> step_keys;
    {
        auto dataset = store.dataset<h5::str, 1>("snapshots/.steps");
        step_keys.resize(dataset.shape().size());
        dataset.read(step_keys.data(), dataset.shape());
    }

    // Start analysis.
    md::periodic_box const box {
        .x_period = box_size,
        .y_period = box_size,
        .z_period = box_size,
    };
    distance_histogram histogram{config.bin_width, config.max_distance, box};

    for (auto const& step_key : step_keys) {
        auto const frame_path = "snapshots/" + step_key;

        // Load and select particle points.
        std::vector<md::point> points;
        {
            auto dataset = store.dataset<float, 2>(frame_path + "/positions");
            points.resize(dataset.shape().dims[0]);
            dataset.read(
                reinterpret_cast<md::scalar*>(points.data()), dataset.shape()
            );
        }

        for (std::size_t i = 0; i < particle_selection.size(); i++) {
            points[i] = points[particle_selection[i]];
        }
        points.resize(particle_selection.size());

        // Calculate RDF. We reuse histogram to reduce memory allocation.
        histogram.clear();
        histogram.update(points);

        for (std::size_t i = 0; i < histogram.size(); i++) {
            auto const density = histogram.density(i);
            auto const posterior = density / expected_density;

            if (i > 0) {
                std::cout << '\t';
            }
            std::cout << posterior;
        }
        std::cout << '\n';
    }
}
