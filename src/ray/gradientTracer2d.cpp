#include <iostream>
#include <limits>
#include <CL/sycl.hpp>
#include <algorithm>
#include <map>
#include <vector>
#ifndef NDEBUG
#include <cassert>
#endif
#if __has_include(<pstl/execution>)
   #include <pstl/execution>
   #include <pstl/algorithm>
   #define USE_PSTL
#endif
//#if __has_include(<oneapi/dpl/execution>)
//   #include <oneapi/dpl/execution>
//   #include <oneapi/dpl/algorithm>
//   #define ONE_API
//#else
//#endif
#include "eikonalxx/ray/gradientTracer2d.hpp"
#include "eikonalxx/ray/path2d.hpp"
#include "eikonalxx/ray/segment2d.hpp"
#include "eikonalxx/ray/point2d.hpp"
#include "eikonalxx/abstractBaseClass/solver2d.hpp"
#include "eikonalxx/geometry2d.hpp"
#include "eikonalxx/source2d.hpp"
#include "eikonalxx/station2d.hpp"
#include "private/grid.hpp"

using namespace EikonalXX::Ray;

namespace
{

template<typename T>
T bilinear(const T x, const T z,
           const T x1, const T x2,
           const T z1, const T z2,
           const T f11, const T f12,
           const T f21, const T f22,
           const T dxi, const T dzi)
{
    T x2mx = x2 - x;
    T xmx1 = x - x1;
    T fxz1 = dxi*(x2mx*f11 + xmx1*f12);
    T fxz2 = dxi*(x2mx*f12 + xmx1*f22);

    T z2mz = z2 - z;
    T zmz1 = z - z1; 
    return dzi*(z2mz*fxz1 + zmz1*fxz2);
}

template<typename T>
T bilinear(const T x, const T z,
           const T x1, const T x2, 
           const T z1, const T z2, 
           const T f11, const T f12,
           const T f21, const T f22)
{
    T dxi = static_cast<T> (1./(x2 - x1));
    T dzi = static_cast<T >(1./(z2 - z1));
    return ::bilinear(x, z, 
                      x1, x2,
                      z1, z2,
                      f11, f12,
                      f21, f22);
}

/// @result The point where this line intersects the next cell wall.
/// The idea for intersection is pretty straightforward.  We have a 
/// point and a gradient, hence we are moving along the linear 
/// trajectory
///   \textbf{f}(t) = \textbf{o} + \textbf{g} t 
/// where \textbf{o} is the origin of the vector, \textbf{g} is 
/// direction given by the gradient field, and t is some unknown value. 
/// Now, we find where \textbf{v} intersects the edge of the cell.  
/// Note, the four (left, right, bottom, top) cell faces are:
///   F_l = x_0
///   F_r = x_1
///   F_b = z_0 
///   F_t = z_1
/// Really, all we are now doing is finding the point of intersection
/// with some existence conditions:
///   Left:   o_x + g_x t_l = x_0 -> t_l = (o_x - x_0)/g_x; g_x < 0
///   Right:  o_x + g_x t_r = x_1 -> t_r = (o_x - x_1)/g_x; g_x > 0
///   Bottom: o_z + g_z t_b = z_0 -> t_b = (o_z - z_0)/g_z; g_z < 0
///   Top:    o_z + g_z t_t = z_1 -> t_t = (o_z - z_1)/g_z; g_z > 0
/// Barring strict equality of a gradient component with zero, we
/// will have two candidate solutions for t; one which intersects the
/// appropriate x plane and one which intersects the appropriate z plane.
/// The t we should pick should be the t that yields the shortest line
/// - i.e., intersects the closest plane.  Therefore, if 
///     | \textbf{o} + \textbf{g} t_i | < | \textbf{o} + \textbf{g} t_j |
/// we would select t_i, otherwise, we select t_j. Following this 
template<typename T>
std::pair<T, T> intersect(const T ox, const T oz,
                          const T gx, const T gz,
                          const T dx, const T dz)
{
    T ti{1.e10};
    T tj{1.e10};
    // Horizontal planes
    T z0 = oz - dz;
    T z1 = oz + dz;
    // Vertical planes
    T x0 = ox - dx;
    T x1 = ox + dx;
    /// Existence
    if (gx > 0) 
    {
        ti = (ox - x1)/gx;
    }
    else if (gx < 0)
    {
        ti = (ox - x0)/gx;
    }
    if (gz > 0)
    {
        tj = (oz - z1)/gz;
    }
    else if (gz > 0)
    {
        tj = (oz - z0)/gz;
    }
    // Make my lines
    T fxi = ox + gx*ti;
    T fzi = oz + gz*ti;
    // Distance squared is good enough
    T li = (fxi - ox)*(fxi - ox) + (fzi - oz)*(fzi - oz);
    T fxj = ox + gx*tj;
    T fzj = oz + gz*tj;
    T lj = (fxj - ox)*(fxj - ox) + (fzj - oz)*(fzj - oz);
    if (li < lj)
    {
        return std::pair {fxi, fzi};
    }
    else
    {
        return std::pair {fxj, fzj};
    }
}

}

class GradientTracer2D::GradientTracer2DImpl
{
public:
    EikonalXX::Geometry2D mGeometry;
    std::vector<EikonalXX::Station2D> mStations;
    std::vector<EikonalXX::Ray::Path2D> mPaths;
    bool mInitialized{false};
};

/// Constructor
GradientTracer2D::GradientTracer2D() :
    pImpl(std::make_unique<GradientTracer2DImpl> ())
{
}

/// Copy constructor
GradientTracer2D::GradientTracer2D(const GradientTracer2D &tracer)
{
    *this = tracer;
}

/// Move constructor
GradientTracer2D::GradientTracer2D(GradientTracer2D &&tracer) noexcept
{
    *this = std::move(tracer);
}

/// Destructor
GradientTracer2D::~GradientTracer2D() = default;

/// Copy assignment
GradientTracer2D& GradientTracer2D::operator=(const GradientTracer2D &tracer)
{
    if (&tracer == this){return *this;}
    pImpl = std::make_unique<GradientTracer2DImpl> (*tracer.pImpl);
    return *this;
}

/// Move assignment
GradientTracer2D&
GradientTracer2D::operator=(GradientTracer2D &&tracer) noexcept
{
    if (&tracer == this){return *this;}
    pImpl = std::move(tracer.pImpl);
    return *this;
}

/// Have stations?
bool GradientTracer2D::haveStations() const noexcept
{
    return !pImpl->mStations.empty();
}

/// Initialized?
bool GradientTracer2D::isInitialized() const noexcept
{
    return pImpl->mInitialized;
}

/// Trace
template<typename T>
void GradientTracer2D::trace(
    const EikonalXX::AbstractBaseClass::ISolver2D<T> &solver)
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    if (!solver.haveVelocityModel())
    {
        throw std::invalid_argument("Velocity model not set on solver");
    }
    if (!solver.haveSource())
    {
        throw std::invalid_argument("Source not on solver");
    }
    if (!solver.haveTravelTimeGradientField())
    {
        throw std::invalid_argument("Travel time gradient field not set");
    }
    // Grid properties
    int nGridX = pImpl->mGeometry.getNumberOfGridPointsInX();
    int nGridZ = pImpl->mGeometry.getNumberOfGridPointsInZ();
    double dx = pImpl->mGeometry.getGridSpacingInX();
    double dz = pImpl->mGeometry.getGridSpacingInZ();
    double x0 = pImpl->mGeometry.getOriginInX();
    double z0 = pImpl->mGeometry.getOriginInZ();
    // Source properties
    auto source = solver.getSource();
    auto iSourceCell = source.getCell();
    auto xs = source.getOffsetInX();
    auto zs = source.getOffsetInZ();
    // Set the source point (we'll always use this when tracing)
    Point2D sourcePoint{x0 + xs, z0 + zs}; 
    const T* gradientPtr = solver.getTravelTimeGradientFieldPointer();
    // Loop on receivers
    for (const auto &station : pImpl->mStations)
    {
        auto iStationCellX = station.getCellInX();
        auto iStationCellZ = station.getCellInZ();
        // Define the station piont (we'll also always use this when tracing)
        auto xr = station.getOffsetInX();
        auto zr = station.getOffsetInZ();
        Point2D stationPoint {x0 + xr, z0 + zr};
        auto iStationCell = station.getCell();
        // Deal with an algorithm breakdown where the source/station are in
        // the same cell
        if (iSourceCell == iStationCell)
        {
            Segment2D segment;
            segment.setStartAndEndPoint(std::pair {sourcePoint, stationPoint});
            segment.setSlowness(
                solver.getSlowness(station.getCellInX(), station.getCellInZ()));
            segment.setVelocityModelCellIndex(iStationCell); 
            continue; 
        }
        // Business as usual - march up the gradient
        bool mConverged{false};
        for (int iSegment = 0;
             iSegment < std::numeric_limits<int>::max(); ++iSegment)
        {
            // 
        }
        if (!mConverged)
        {
            std::cerr << "Ray did not converge to source" << std::endl;
        }
    }
}

///--------------------------------------------------------------------------///
///                            Template Instantiation                        ///
///--------------------------------------------------------------------------///
template void EikonalXX::Ray::GradientTracer2D::trace(
    const EikonalXX::AbstractBaseClass::ISolver2D<double> &solver);
template void EikonalXX::Ray::GradientTracer2D::trace(
    const EikonalXX::AbstractBaseClass::ISolver2D<float> &solver);
