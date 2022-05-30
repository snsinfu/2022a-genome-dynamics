#pragma once

// This header file defines `h5::buffer_traits` specializations for saving
// simulation data to an HDF5 file.

#include <cstddef>

#include <h5.hpp>
#include <md.hpp>

#include "loops/loop_simulator.hpp"


template<>
struct h5::buffer_traits<md::array_view<md::point const> const>
{
    using buffer_type = md::array_view<md::point const>;
    using value_type = md::scalar;
    static constexpr std::size_t rank = 2;
    static constexpr std::size_t dimension = 3;

    static h5::shape<rank> shape(buffer_type const& buffer)
    {
        return {buffer.size(), dimension};
    }

    static value_type const* data(buffer_type const& buffer)
    {
        return &buffer.data()->x;
    }
};


template<>
struct h5::buffer_traits<md::array_view<loop_pair const> const>
{
    using buffer_type = md::array_view<loop_pair const>;
    using value_type = std::size_t;
    static constexpr std::size_t rank = 2;
    static constexpr std::size_t dimension = 3;

    static h5::shape<rank> shape(buffer_type const& buffer)
    {
        return {buffer.size(), dimension};
    }

    static value_type const* data(buffer_type const& buffer)
    {
        return &buffer.data()->start;
    }
};


template<>
struct h5::buffer_traits<md::array_view<int const> const>
{
    using value_type = int;
    using buffer_type = md::array_view<value_type const>;
    static constexpr std::size_t rank = 1;

    static h5::shape<rank> shape(buffer_type const& buffer)
    {
        return {buffer.size()};
    }

    static value_type const* data(buffer_type const& buffer)
    {
        return buffer.data();
    }
};
