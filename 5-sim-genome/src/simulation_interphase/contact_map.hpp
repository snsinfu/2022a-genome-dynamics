#pragma once

// This module defines the interface of contact_map class. The class computes
// time-integrated contact map of moving points.

#include <cstddef>
#include <cstdint>
#include <functional>
#include <vector>

#include <Eigen/Sparse>

#include <md.hpp>


// Class: contact_map
//
// Accumulates time-integrated contact map of moving points.
//
class contact_map
{
    using Count = std::uint32_t;
    using Matrix = Eigen::SparseMatrix<Count, Eigen::RowMajor, std::int32_t>;
    using Triplet = Eigen::Triplet<Count, std::int32_t>;

public:
    // Function: set_contact_distance
    //
    // Sets contact distance used in following update.
    //
    void set_contact_distance(md::scalar dist);

    // Function: contact_distance
    //
    // Returns the current contact distance.
    //
    md::scalar contact_distance() const;

    // Function: clear
    //
    // Clears contact map in-place.
    //
    void clear();

    // function: update
    //
    // Computes contact map of given points and adds to the ensemble.
    //
    void update(md::array_view<const md::point> points);

    // Function: accumulate
    //
    // Returns (i,j,v)-style contact map where i and j are indices and v is the
    // number of contacts.
    //
    std::vector<std::array<Count, 3>> accumulate() const;

private:
    md::scalar _contact_distance = 0;
    Matrix     _contact_matrix;

    std::vector<Triplet> _triplet_buffer;
    Matrix               _matrix_buffer;
};
