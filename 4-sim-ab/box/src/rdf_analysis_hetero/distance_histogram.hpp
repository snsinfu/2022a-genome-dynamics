#pragma once

#include <cstddef>
#include <vector>

#include <md.hpp>


class distance_histogram
{
public:
    explicit distance_histogram(
        md::scalar bin_width,
        md::scalar max_distance,
        md::periodic_box const& box
    );

    void update(
        md::array_view<md::point const> center_points,
        md::array_view<md::point const> target_points
    );
    void        clear();
    std::size_t size() const;
    double      frequency(std::size_t i) const;
    double      density(std::size_t i) const;

private:
    md::scalar                              _bin_width;
    md::scalar                              _max_distance;
    md::periodic_box                        _box;
    md::neighbor_searcher<md::periodic_box> _searcher;
    std::vector<double>                     _bin_freqs;
    std::vector<double>                     _bin_volumes;
};
