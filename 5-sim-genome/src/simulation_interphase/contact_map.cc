#include <array>
#include <cassert>
#include <cstdint>
#include <utility>
#include <vector>

#include <Eigen/Sparse>
#include <md.hpp>

#include "contact_map.hpp"


void contact_map::set_contact_distance(md::scalar dist)
{
    _contact_distance = dist;
}


md::scalar contact_map::contact_distance() const
{
    return _contact_distance;
}


void contact_map::clear()
{
    _contact_matrix.setZero();
}


void contact_map::update(md::array_view<const md::point> points)
{
    auto const matrix_dim = static_cast<Matrix::Index>(points.size());
    if (_contact_matrix.rows() < matrix_dim) {
        _contact_matrix.resize(matrix_dim, matrix_dim);
    }

    struct triplet_output_iterator
    {
        std::vector<Triplet>& output;

        triplet_output_iterator operator++(int)
        {
            return *this;
        }

        triplet_output_iterator& operator*()
        {
            return *this;
        }

        void operator=(std::pair<md::index, md::index> const& pair)
        {
            output.push_back({
                static_cast<std::int32_t>(pair.first),
                static_cast<std::int32_t>(pair.second),
                1
            });
        }
    };

    _triplet_buffer.clear();

    md::neighbor_searcher<md::open_box> searcher{{}, _contact_distance};
    searcher.set_points(points);
    searcher.search(triplet_output_iterator{_triplet_buffer});

    _matrix_buffer.resize(_contact_matrix.rows(), _contact_matrix.cols());
    _matrix_buffer.setFromTriplets(_triplet_buffer.begin(), _triplet_buffer.end());

    _contact_matrix += _matrix_buffer;
}


std::vector<std::array<std::uint32_t, 3>> contact_map::accumulate() const
{
    std::vector<std::array<std::uint32_t, 3>> contacts;
    contacts.reserve(static_cast<std::size_t>(_contact_matrix.nonZeros()));

    for (Matrix::Index i = 0; i < _contact_matrix.outerSize(); ++i) {
        for (Matrix::InnerIterator it{_contact_matrix, i}; it; ++it) {
            contacts.push_back({
                static_cast<std::uint32_t>(it.row()),
                static_cast<std::uint32_t>(it.col()),
                it.value()
            });
        }
    }

    return contacts;
}
