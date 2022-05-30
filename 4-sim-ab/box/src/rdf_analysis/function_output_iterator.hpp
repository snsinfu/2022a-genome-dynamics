#pragma once


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

    function_output_iterator& operator++()
    {
        return *this;
    }

    function_output_iterator operator++(int)
    {
        return *this;
    }

    template<typename T>
    void operator=(T const& value)
    {
        _func(value);
    }

private:
    F _func;
};

template<typename F>
function_output_iterator<F> make_function_output_iterator(F const& func)
{
    return function_output_iterator<F>{func};
}
