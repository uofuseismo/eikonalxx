#ifndef PRIVATE_PAD_HPP
#define PRIVATE_PAD_HPP
namespace
{
/// @brief Computes the matrix padding length so that we get bit aligned memory.
/// @param[in] n              The array length.  This must be positive.
/// @param[in] precisionSize  The size of the precision, ex: sizeof(float) = 4.
/// @param[in] alignment      The byte alignment.  This should be a power of 2.
/// @result The padded array length.  This will be greater than or equal to n.
inline int padLength(const int n,
                     const size_t precisionSize = sizeof(double),
                     const int alignment = 64)
{
    auto size = static_cast<int> (precisionSize);
    int padLength = 0;
    auto xmod = (n*size)%alignment;
    if (xmod != 0){padLength = (alignment - xmod)/size;}
    auto nptsPadded = n + padLength;
    return nptsPadded;
}
}
#endif
