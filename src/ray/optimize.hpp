#ifndef PRIVATE_OPTIMIZE_HPP
#define PRIVATE_OPTIMIZE_HPP
/// @brief Utilities for optimizing ray paths in ray shooting.
#include <cmath>
#include <vector>
#include <functional>
#ifndef NDEBUG
#include <cassert>
#endif
#include "layerTracer.hpp"
#include "eikonalxx/ray/path2d.hpp"
namespace
{
/// @result 
void computeTakeOffAngles(
    const int n, const double theta0, const double theta1,
    std::vector<double> *takeOffAngles)
{
#ifndef NDEBUG
    assert(n > 0);
#endif
    if (static_cast<int> (takeOffAngles->size()) != n)
    {
        takeOffAngles->resize(n);
    }
    if (n == 1)
    {
        takeOffAngles->at(0) = 0.5*(theta0 + theta1);
    }
    auto dTheta = (theta1 - theta0)/(n - 1);
    auto *__restrict__ takeOffAnglePointer = takeOffAngles->data();
    for (int i = 0; i < n; ++i)
    {
        takeOffAnglePointer[i] = theta0 + dTheta*i;
    }
}

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

/*
void optimizeUpGoing()
{

}
*/

std::vector<EikonalXX::Ray::Path2D>
   optimizeCriticallyRefracted(
    const std::vector<double> &augmentedInterfaces,
    const std::vector<double> &augmentedSlownesses,
    const int sourceLayer,
    const double sourceDepth,
    const int stationLayer,
    const double stationDepth,
    const double stationOffset)
{
    constexpr double rayHitTolerance{1}; // 1 meter
    std::vector<EikonalXX::Ray::Path2D> rayPaths;
    // Flip the problem around
    if (stationDepth > sourceDepth)
    {
        std::cout << "Station below source for critical refraction" << std::endl;
        rayPaths = ::optimizeCriticallyRefracted(augmentedInterfaces,
                                                 augmentedSlownesses,
                                                 stationLayer,
                                                 sourceDepth,
                                                 stationLayer,
                                                 stationDepth,
                                                 stationOffset);
        for (auto &rayPath : rayPaths)
        {
            rayPath.reverse();
        }
        return rayPaths;
    }
    // The idea here is pretty simple.  Basically, the trick is to spend
    // the least amount of time possible in a slower layer.  To do that
    // we use the smallest angle possible - which is the critical angle.
    // Then we just use some algebra to relate that angle for each layer
    // to the other layers.  Since:
    //    sin(i_c) = V_{i}/V_{i+1} -> i_c = asin(V_{i}/V_{i+1})
    // Then what's this angle in the previous layer?
    //    sin(i)/V_{i-1} = sin(i_c)/V_i
    //                   = (V_i/V_{i+1})/V_i
    //                   = 1/V{i+1}
    //    sin(i)/V_{i-1} = 1/V_{i+1}
    // -> i = asin(V_{i-1}/V_{i+1})
    // And following the recursion down we get, in general, for a velocity
    // model where the velocities increase with depth:
    // -> i = asin(V_{sourceLayer}/V_{i+1})
    // Step 1: Tabulate all the critical angles
    std::vector<double> takeOffAngles;
    int nLayers = static_cast<int> (augmentedSlownesses.size()) - 1;
    auto lastLayer = nLayers - 1;
    takeOffAngles.reserve(nLayers);
    for (int layer = sourceLayer + 1; layer < nLayers; ++layer)
    {
        if (augmentedSlownesses[sourceLayer] < augmentedSlownesses[layer])
        {
            std::cerr << "Velocity inversion exists!" << std::endl;
        }
        auto angle = ::computeCriticalAngle(augmentedSlownesses[sourceLayer],
                                            augmentedSlownesses[layer]);
        // Pad just a touch in case we get burned by inexact math
        takeOffAngles.push_back(angle*(180/M_PI) + 1.e-10);
    }
    // Shoot these rays
    for (const auto &takeOffAngle : takeOffAngles)
    {
        constexpr bool allowCriticalRefractions{true};
        constexpr bool keepOnlyHits{true}; 
        auto paths = ::shoot(takeOffAngle,
                             augmentedInterfaces,
                             augmentedSlownesses,
                             sourceLayer,
                             sourceDepth,
                             stationLayer,
                             stationOffset,
                             stationDepth,
                             allowCriticalRefractions,
                             rayHitTolerance,
                             keepOnlyHits,
                             lastLayer);
        if (!paths.empty())
        {
            rayPaths.insert(rayPaths.end(), paths.begin(), paths.end());
        }
    }
    // Sort them
    if (!rayPaths.empty())
    {
        std::sort(rayPaths.begin(), rayPaths.end(),
                  [](const EikonalXX::Ray::Path2D &lhs,
                     const EikonalXX::Ray::Path2D &rhs)
                 {
                     return lhs.getTravelTime() < rhs.getTravelTime();
                 });
        //std::cout << rayPaths[0].getTravelTime() << std::endl;
    }
    return rayPaths;
}
    

/// Slownesses - augmented slownesses
void optimizeFirstArrivingDownGoing(
    const std::vector<double> &augmentedInterfaces,
    const std::vector<double> &augmentedSlownesses,
    const double sourceOffset,
    const int sourceLayer,
    const double sourceDepth,
    const int stationLayer,
    const double stationDepth,
    const double stationOffset,
    const double rayHitTolerance)
{
    constexpr bool allowCriticalRefractions{false};
    auto nLayers = static_cast<int> (augmentedSlownesses.size()) - 2;
    if (sourceLayer >= nLayers){return;}
    // Compute an upper bound on angles
    double maximumTakeOffAngle = M_PI_2 - 1.e-4;
    for (int layer = sourceLayer + 1; layer < nLayers; ++layer)
    {
        auto criticalAngle
            = ::computeCriticalAngle(augmentedSlownesses[sourceLayer],
                                     augmentedSlownesses[layer]);
        maximumTakeOffAngle = std::min(maximumTakeOffAngle, criticalAngle);
    }
    constexpr double minimumTakeOffAngle{0};
    maximumTakeOffAngle = std::max(1.e-4, maximumTakeOffAngle - 1.e-8);
    int nAngles = static_cast<int> (maximumTakeOffAngle) + 1;
    std::vector<double> takeOffAngles;
    ::computeTakeOffAngles(nAngles, minimumTakeOffAngle, maximumTakeOffAngle,
                           &takeOffAngles);
    // Compute the first batch of ray paths in the coarse grid search
    auto lastLayer = static_cast<int> (augmentedInterfaces.size()) - 2;
    for (int i = static_cast<int> (takeOffAngles.size()) - 1; i >= 0; --i)
    {
        constexpr bool keepOnlyHits{false};
        auto paths = ::shoot(takeOffAngles[i],
                             augmentedInterfaces,
                             augmentedSlownesses,
                             sourceLayer,
                             sourceDepth,
                             stationLayer,
                             stationOffset,
                             stationDepth,
                             allowCriticalRefractions,
                             rayHitTolerance,
                             keepOnlyHits,
                             lastLayer);
        if (!paths.empty())
        {
            
        }
    }
    // If the offset exceeds this then quit
    
}

}
#endif
