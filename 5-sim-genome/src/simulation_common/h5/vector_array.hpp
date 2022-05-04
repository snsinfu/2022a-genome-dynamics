#pragma once

// This header defines a HighFive trait for storing std::vector<std::array> as
// two-dimensional CArray.

#include <array>
#include <cstddef>
#include <string>
#include <vector>
#include <cassert>

#include <highfive/H5Attribute.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5Group.hpp>
#include <highfive/H5PropertyList.hpp>


namespace H5 = HighFive;


namespace HighFive::details
{
    // Adapt std::vector<std::array<T, D>> as (*,D)-CArray.
    template<typename T, std::size_t D>
    struct data_converter<std::vector<std::array<T, D>>, void>
    {
        using row_type = std::array<T, D>;
        using array_type = std::vector<row_type>;

        static_assert(
            sizeof(row_type) == sizeof(T[D]),
            "This code assumes no padding in std::array"
        );

        data_converter(array_type& array, DataSpace& space, std::size_t dim = 0)
        {
            const auto dims = space.getDimensions();
            assert(dims.size() == 2);
            assert(dims[1] == D);
            (void) array;
            (void) dim;
            _rows = dims[0];
        }

        T* transform_read(array_type& array)
        {
            array.resize(_rows);
            return array[0].data();
        }

        T* transform_write(array_type& array)
        {
            return array[0].data();
        }

        void process_result(array_type& array)
        {
            // Nothing to do.
            (void) array;
        }

    private:
        std::size_t _rows;
    };
}
