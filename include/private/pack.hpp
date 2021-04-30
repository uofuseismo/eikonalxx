#ifndef PRIVATE_PACK_HPP
#define PRIVATE_PACK_HPP
#include <vector>
#include <algorithm>
#include <string>
#include <bit>

namespace
{

std::string fillBlanksWithUnderscores(const std::string &s) 
{
    auto result = s;
    std::replace_if(result.begin(), result.end(), ::isspace, '_');
    return result;
}

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

void pack(const size_t n, const float *__restrict__ x,
          char *__restrict__ c4) 
{
    for (size_t i = 0; i < n; ++i)
    {   
        pack(x[i], &c4[4*i]);
    }   
}

std::vector<char> pack(const std::vector<float> &x) 
{
    std::vector<char> cx(x.size()*4);
    pack(x.size(), x.data(), cx.data());
    return cx; 
}

}
#endif
