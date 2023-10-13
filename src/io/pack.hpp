#ifndef PRIVATE_PACK_HPP
#define PRIVATE_PACK_HPP
#include <vector>
#include <algorithm>
#include <string>

namespace
{

std::string fillBlanksWithUnderscores(const std::string &s) 
{
    auto result = s;
    std::replace_if(result.begin(), result.end(), ::isspace, '_');
    return result;
}

[[maybe_unused]]
void pack(const float f4, char c4[4])
{
    auto cTemp = reinterpret_cast<const char *> (&f4);
    if constexpr (std::endian::native == std::endian::little)
    {
        c4[0] = cTemp[3];
        c4[1] = cTemp[2];
        c4[2] = cTemp[1];
        c4[3] = cTemp[0];
    }
    else
    {   
        std::copy(cTemp, cTemp + 4, c4);
    }   
}

[[maybe_unused]]
void pack(const double f8, char c4[4])
{
    return ::pack(static_cast<float> (f8), c4);
}

[[maybe_unused]]
void pack(const size_t n, const float *__restrict__ x,
          char *__restrict__ c4) 
{
    for (size_t i = 0; i < n; ++i)
    {   
        ::pack(x[i], &c4[4*i]);
    }   
}

[[maybe_unused]]
void pack(const size_t n, const double *__restrict__ x,
          char *__restrict__ c4) 
{
    for (size_t i = 0; i < n; ++i)
    {
        ::pack(x[i], &c4[4*i]);
    }
}

[[maybe_unused]]
std::vector<char> pack(const std::vector<float> &x) 
{
    std::vector<char> cx(x.size()*4);
    ::pack(x.size(), x.data(), cx.data());
    return cx; 
}

[[maybe_unused]]
std::vector<char> pack(const std::vector<double> &x) 
{
    std::vector<char> cx(x.size()*8);
    ::pack(x.size(), x.data(), cx.data());
    return cx; 
}


}
#endif
