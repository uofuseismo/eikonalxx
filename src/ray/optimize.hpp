#ifndef PRIVATE_OPTIMIZE_HPP
#define PRIVATE_OPTIMIZE_HPP
/// @brief Utilities for optimizing ray paths in ray shooting.
#include <cmath>
#include <vector>
#include <functional>
#ifndef NDEBUG
#include <cassert>
#endif
#include "eikonalxx/ray/path2d.hpp"
namespace
{
[[nodiscard]]
/// @param[in] shootRay  A function, that given a take-off angle, shoots a
///                      ray.
/// @param[in] fa  The left-bracketing take-off angle and travel time.
/// @param[in] fb  An intermediate take-off angle and travel time.
/// @param[in] fc  The right-bracketing take-off angle and travel time.
/// @param[in] nIterations  The number of iterations in the optimization.
/// @result The optimal take-off angle and ray path.
std::pair<double, EikonalXX::Ray::Path2D>
    quadraticFitSearch(std::function< std::vector<EikonalXX::Ray::Path2D>
                       (const double)> &shootRay,
                       const std::pair<double, double> &fa,
                       const std::pair<double, double> &fb,
                       const std::pair<double, double> &fc,
                       const int nIterations = 5)
{
    constexpr double third{1./3.};
    auto a  = fa.first;
    auto ya = fa.second;
    auto b  = fb.first;
    auto yb = fb.second;
    auto c  = fc.first;
    auto yc = fc.second;
    auto xStar = third*(a + b + c);
    auto fStar = std::max(ya, std::max(yb, yc));
    EikonalXX::Ray::Path2D fPath;
    for (int k = 0; k < nIterations; ++k)
    {
#ifndef NDEBUG
        assert(a < b);
        assert(b < c);
#endif
        auto a2 = a*a;
        auto b2 = b*b;
        auto c2 = c*c;
        auto xNum = ya*(b2 - c2) + yb*(c2 - a2) + yc*(a2 - b2);
        auto xDen = 2*(ya*(b - c) + yb*(c - a) +  yc*(a - b));
        xStar = third*(a + b + c);
        if (std::abs(xDen) > 1.e-14)
        {
            xStar = xNum/xDen;
        }
        else
        {
            std::cerr << "Division by 0: " << xNum << ","
                      << xDen << std::endl;
        }
        std::vector<EikonalXX::Ray::Path2D> paths;
        try
        {
            paths = shootRay(xStar);
            if (paths.empty())
            {
                throw std::runtime_error("No paths converged an angle: "
                                       + std::to_string(xStar));
            }
        }
        catch (const std::exception &e)
        {
            throw std::runtime_error("Optimization failed: "
                                   + std::string{e.what()});
        }
        fPath = paths.at(0);
        fStar = fPath.getTravelTime();
        if (xStar > b)
        {
            if (fStar > yb) 
            {
                c = xStar;
                yc = fStar; 
            }
            else
            {
                a  = b;
                ya = yb; 
                b  = xStar;
                yb = fStar;
            }   
        }
        else if (xStar < b)
        {
            if (fStar > yb) 
            {
                a = xStar;
                ya = fStar;
            }
            else
            {
                c  = b;
                yc = yb; 
                b  = xStar;
                yb = fStar;
            }
        }
        else
        {
            return std::pair {b, fPath};
        }
    }
    return std::pair {xStar, fPath};
}

}
#endif
