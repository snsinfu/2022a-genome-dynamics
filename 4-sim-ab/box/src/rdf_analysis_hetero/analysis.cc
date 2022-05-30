#include <cmath>
#include <cstddef>
#include <iostream>
#include <stdexcept>
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

    // Select center particles and target particles for RDF calculation. RDF
    // then measures the density of target particle around center particle.
    md::scalar center_a = -1;
    if (config.center_type == "A") {
        center_a = 1;
    }
    if (config.center_type == "B") {
        center_a = 0;
    }
    if (center_a == -1) {
        throw std::runtime_error("invalid center type: '" + config.center_type + "'");
    }

    std::vector<md::index> center_indices;
    std::vector<md::index> target_indices;

    for (md::index i = 0; i < ab_factors.size(); i++) {
        if (std::fabs(ab_factors[i].first - center_a) < 1e-6) {
            center_indices.push_back(i);
        } else {
            target_indices.push_back(i);
        }
    }

    // Calculate expected density of target points. RDF is normalized to this
    // density.
    md::scalar box_size;
    {
        std::string config_json;
        store.dataset<h5::str>("metadata/config").read(config_json);
        auto const simulation_config = nlohmann::json::parse(config_json);
        box_size = simulation_config["box_size"];
    }
    auto const volume = md::power<3>(box_size);
    auto const expected_density = md::scalar(target_indices.size()) / volume;

    // Get snapshot keys to analyze. We load all.
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

        // Load and select center and target points.
        std::vector<md::point> points;
        {
            auto dataset = store.dataset<float, 2>(frame_path + "/positions");
            points.resize(dataset.shape().dims[0]);
            dataset.read(
                reinterpret_cast<md::scalar*>(points.data()), dataset.shape()
            );
        }

        std::vector<md::point> center_points;
        std::vector<md::point> target_points;

        for (auto const i : center_indices) {
            center_points.push_back(points[i]);
        }

        for (auto const i : target_indices) {
            target_points.push_back(points[i]);
        }

        // Calculate RDF. This is the bottleneck. We reuse histogram to reduce
        // memory allocation.
        histogram.clear();
        histogram.update(center_points, target_points);

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
