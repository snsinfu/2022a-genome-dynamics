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


void distance_histogram::update(md::array_view<md::point const> points)
{
    // Neighbor searcher scans unique pairs, so (j,i) is not covered when (i,j)
    // is. Double the weight to count both.
    auto const unit_weight = 2 / md::scalar(points.size());

    _searcher.set_points(points);
    _searcher.search(
        make_function_output_iterator([&](auto pair) {
            auto const [ i, j ] = pair;
            auto const disp = _box.shortest_displacement(points[i], points[j]);
            auto const distance = disp.norm();
            auto const bin_index = std::size_t(distance * (1 / _bin_width));
            if (bin_index >= _bin_freqs.size()) {
                return; // Cut off.
            }
            _bin_freqs[bin_index] += unit_weight;
        })
    );
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
