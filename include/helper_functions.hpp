#ifndef GET_SIZE_HPP
#define GET_SIZE_HPP
#include <iostream>

constexpr inline size_t get_num_elements(float low, float high, float step)
{
    float range = (high - low);
    float N     = range/step;
    return (size_t)N  + 1;
}  

template<typename T>
concept has_size = requires(T matrix)
{
    {matrix.size()} -> std::convertible_to<size_t>;
};
size_t get_size(const has_size auto& matrix)
{
    auto size = static_cast<float>(matrix.size()*sizeof(std::complex<float>));
    return static_cast<size_t>(size);
}

#endif
