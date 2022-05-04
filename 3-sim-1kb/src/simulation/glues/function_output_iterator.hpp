#pragma once

#include <utility>


template<typename F>
class function_output_iterator
{
public:
    explicit function_output_iterator(F const& func)
        : _func{func}
    {
    }

    function_output_iterator& operator*()
    {
        return *this;
    }

    template<typename T>
    void operator=(T&& value)
    {
        _func(std::forward<T>(value));
    }

    function_output_iterator& operator++()
    {
        return *this;
    }

    function_output_iterator& operator++(int)
    {
        return *this;
    }

private:
    F _func;
};