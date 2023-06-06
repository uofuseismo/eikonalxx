#ifndef BILINEAR_HPP
#define BILINEAR_HPP
namespace
{
template<typename T>
T bilinear(const T x, const T z,
           const T x1, const T x2, 
           const T z1, const T z2, 
           const T f00, const T f01,
           const T f10, const T f11,
           const T dxi, const T dzi)
{
    T x2mx = x2 - x;
    T xmx1 = x - x1; 
    T fxz1 = dxi*(x2mx*f00 + xmx1*f10);
    T fxz2 = dxi*(x2mx*f01 + xmx1*f11);

    T z2mz = z2 - z;
    T zmz1 = z - z1; 
    return dzi*(z2mz*fxz1 + zmz1*fxz2);
}
}
#endif
