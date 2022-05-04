#pragma once

#include <md.hpp>


struct particle_data
{
    md::scalar a_factor = 0;
    md::scalar b_factor = 0;
};


inline md::attribute_key<particle_data> particle_data_attribute;
