#include <cmath>
#include <cstddef>

#include <md.hpp>

#include "distance_histogram.hpp"
#include "function_output_iterator.hpp"


namespace
{
    constexpr md::scalar PI = 3.1416;
}


distance_histogram::distance_histogram(
    md::scalar bin_width,
    md::scalar max_distance,
    md::periodic_box const& box
)
    : _bin_width{bin_width}
    , _max_distance{max_distance}
    , _box{box}
    , _searcher{box, max_distance}
    , _bin_freqs(std::size_t(std::ceil(_max_distance / _bin_width)))
{
    for (std::size_t i = 0; i < _bin_freqs.size(); i++) {
        auto r_min = _bin_width * md::scalar(i);
        auto r_max = _bin_width * md::scalar(i + 1);
        if (r_max > max_distance) {
            r_max = max_distance;
        }
        auto const dr3 = md::power<3>(r_max) - md::power<3>(r_min);
        auto const volume = 4 * PI / 3 * dr3;
        _bin_volumes.push_back(volume);
    }
}


void distance_histogram::update(
    md::array_view<md::point const> center_points,
    md::array_view<md::point const> target_points
)
{
    // We want to compute mean histogram. So, each center point contributes to
    // frequency by 1/N where N = center_points.size().
    auto const unit_weight = 1 / md::scalar(center_points.size());

    _searcher.set_points(target_points);

    for (auto const& center : center_points) {
        _searcher.query(
            center,
            make_function_output_iterator([&](auto i) {
                auto const disp = _box.shortest_displacement(center, target_points[i]);
                auto const distance = disp.norm();
                auto const bin_index = std::size_t(distance * (1 / _bin_width));
                if (bin_index >= _bin_freqs.size()) {
                    return; // Cut off.
                }
                _bin_freqs[bin_index] += unit_weight;
            })
        );
    }
}


void distance_histogram::clear()
{
    for (auto& freq : _bin_freqs) {
        freq = 0;
    }
}


std::size_t distance_histogram::size() const
{
    return _bin_freqs.size();
}


double distance_histogram::frequency(std::size_t i) const
{
    return _bin_freqs[i];
}

double distance_histogram::density(std::size_t i) const
{
    return _bin_freqs[i] / _bin_volumes[i];
}
